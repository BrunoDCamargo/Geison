import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.config import ConservationConfig
from qpcr_pipeline.conservation import ReferenceAnnotation, WindowConservation
from qpcr_pipeline.report_html import render_conservation_html


DATA_MARKER = '<script id="geison-report-data" type="application/json">'


def window(
    start,
    mean_conservation,
    minimum_conservation,
    mean_coverage,
    mean_entropy_bits,
):
    return WindowConservation(
        reference_start=start,
        reference_end=start + 49,
        position_count=50,
        mean_conservation=mean_conservation,
        minimum_conservation=minimum_conservation,
        mean_coverage=mean_coverage,
        mean_gap_frequency=1.0 - mean_coverage,
        mean_entropy_bits=mean_entropy_bits,
    )


def embedded_data(html):
    payload_start = html.index(DATA_MARKER) + len(DATA_MARKER)
    payload_end = html.index("</script>", payload_start)
    return json.loads(html[payload_start:payload_end])


def embedded_data_source(html):
    payload_start = html.index(DATA_MARKER) + len(DATA_MARKER)
    payload_end = html.index("</script>", payload_start)
    return html[payload_start:payload_end]


def report_javascript(html):
    payload_end = html.index("</script>", html.index(DATA_MARKER) + len(DATA_MARKER))
    script_start = html.index("<script>", payload_end) + len("<script>")
    script_end = html.index("</script>", script_start)
    return html[script_start:script_end]


def javascript_function(source, name):
    start = source.index(f"function {name}(")
    opening_brace = source.index("{", start)
    depth = 0
    for index in range(opening_brace, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[start : index + 1]
    raise AssertionError(f"Unterminated JavaScript function: {name}")


class ConservationReportHtmlTests(unittest.TestCase):
    def render(self, *, target_name="target", windows=(), annotations=()):
        return render_conservation_html(
            target_name=target_name,
            reference_id="reference-1",
            sequence_count=17,
            config=ConservationConfig(enabled=True, window_size=50, step_size=10),
            windows=tuple(windows),
            annotations=tuple(annotations),
        )

    def test_report_is_self_contained_interactive_and_safely_embeds_data(self):
        unsafe_target = "target </script><img src=x onerror=alert(1)> & \u2028\u2029"
        unsafe_label = "gene <b>alpha</b> </script> & \u2028\u2029"
        windows = (
            window(1, 0.95, 0.85, 1.0, 0.2),
            window(11, 0.90, 0.75, 0.9, 0.4),
        )
        annotations = (ReferenceAnnotation("gene", 2, 35, 1, unsafe_label),)

        html = self.render(
            target_name=unsafe_target, windows=windows, annotations=annotations
        )
        data = embedded_data(html)
        payload_start = html.index(DATA_MARKER) + len(DATA_MARKER)
        payload_end = html.index("</script>", payload_start)
        embedded_source = html[payload_start:payload_end]

        self.assertEqual(data["identity"]["targetName"], unsafe_target)
        self.assertEqual(data["annotations"][0][4], unsafe_label)
        self.assertNotIn("</script", embedded_source.lower())
        self.assertNotIn("</script><img", html.lower())
        self.assertNotIn("https://", html)
        self.assertNotIn("http://", html)
        self.assertNotIn('src="//', html)
        self.assertIn("<canvas", html)
        self.assertIn('id="reset-zoom"', html)
        self.assertIn('id="hover-details"', html)
        self.assertIn('id="top-windows"', html)
        self.assertIn('addEventListener("wheel"', html)
        self.assertIn('addEventListener("pointermove"', html)
        self.assertIn('addEventListener("pointerdown"', html)
        self.assertIn("textContent", html)

    def test_report_declares_an_inline_favicon_without_a_network_request(self):
        html = self.render(windows=(window(1, 0.95, 0.85, 1.0, 0.2),))

        self.assertIn('<link rel="icon" href="data:,">', html)
        self.assertNotIn('href="/favicon.ico"', html)

    def test_template_marker_text_round_trips_without_post_serialization_replacement(self):
        marker_text = "target __EMPTY_HIDDEN__ annotation"

        html = self.render(
            target_name=marker_text,
            windows=(window(1, 0.95, 0.85, 1.0, 0.2),),
            annotations=(ReferenceAnnotation("gene", 1, 25, 1, marker_text),),
        )
        data = embedded_data(html)

        self.assertEqual(data["identity"]["targetName"], marker_text)
        self.assertEqual(data["annotations"][0][4], marker_text)

    def test_template_marker_text_round_trips_when_empty_state_is_visible(self):
        marker_text = "empty __EMPTY_HIDDEN__ report"

        html = self.render(
            target_name=marker_text,
            windows=(),
            annotations=(ReferenceAnnotation("gene", 1, 25, 1, marker_text),),
        )
        data = embedded_data(html)

        self.assertEqual(data["identity"]["targetName"], marker_text)
        self.assertEqual(data["annotations"][0][4], marker_text)
        self.assertIn('id="empty-state" role="status" aria-live="polite">', html)

    def test_peak_and_trace_drawing_is_clipped_before_axes_are_drawn(self):
        html = self.render(windows=(window(1, 0.95, 0.85, 1.0, 0.2),))

        save = html.index("context.save();", html.index("function draw()"))
        clip = html.index("context.clip();", save)
        peak = html.index("context.fillRect(xFor(item[0])", clip)
        conservation = html.index('drawTrace(items,3,"#1769aa")', peak)
        coverage = html.index('drawTrace(items,5,"#d65f00")', conservation)
        restore = html.index("context.restore();", coverage)
        axes = html.index('context.fillText(formatCoordinate(viewStart)', restore)

        self.assertLess(save, clip)
        self.assertLess(clip, peak)
        self.assertLess(peak, conservation)
        self.assertLess(conservation, coverage)
        self.assertLess(coverage, restore)
        self.assertLess(restore, axes)

    def test_top_windows_are_ranked_by_the_scientific_tie_breakers(self):
        windows = (
            window(1, 0.80, 0.70, 0.90, 0.50),
            window(2, 0.95, 0.80, 0.80, 0.40),
            window(3, 0.95, 0.80, 0.90, 0.60),
            window(4, 0.95, 0.80, 0.90, 0.20),
            window(5, 0.95, 0.80, 0.90, 0.20),
            window(6, 0.95, 0.70, 1.00, 0.10),
            window(7, 0.90, 0.60, 1.00, 0.10),
            window(8, 0.85, 0.80, 1.00, 0.10),
            window(9, 0.75, 0.70, 1.00, 0.10),
            window(10, 0.70, 0.70, 1.00, 0.10),
            window(11, 0.65, 0.60, 1.00, 0.10),
            window(12, 0.60, 0.60, 1.00, 0.10),
        )

        data = embedded_data(self.render(windows=windows))

        self.assertEqual(
            [item[0] for item in data["topWindows"]],
            [4, 5, 3, 2, 6, 7, 8, 1, 9, 10],
        )

    def test_report_exposes_summary_legend_instructions_and_hover_fields(self):
        html = self.render(
            target_name="Respiratory target",
            windows=(window(1, 0.95, 0.85, 1.0, 0.2),),
            annotations=(ReferenceAnnotation("CDS", 5, 25, -1, "capsid"),),
        )

        for visible_text in (
            "Target",
            "Reference",
            "Discovery sequences",
            "Windows",
            "Conservation",
            "Coverage",
            "Zoom",
            "Pan",
            "Reset",
            "Reference annotations",
            "Interval",
            "Mean conservation",
            "Minimum conservation",
            "Mean coverage",
            "Mean gaps",
            "Mean entropy",
            "Overlapping annotations",
        ):
            with self.subTest(visible_text=visible_text):
                self.assertIn(visible_text, html)
        self.assertIn('row.addEventListener("click"', html)
        self.assertIn("window.devicePixelRatio", html)
        self.assertNotIn("<svg", html.lower())

    def test_hover_runtime_searches_by_window_centers(self):
        node = shutil.which("node")
        if node is None:
            self.skipTest("Node.js is unavailable for report JavaScript verification")
        windows = tuple(
            window(start, 0.95, 0.85, 1.0, 0.2)
            for start in range(1, 202, 10)
        )
        html = self.render(windows=windows)
        source = report_javascript(html)
        nearest_window = javascript_function(source, "nearestWindow")
        runtime_source = (
            'const windows=JSON.parse(require("fs").readFileSync(0,"utf8"));'
            "const viewStart=1;const viewEnd=250;"
            f"{nearest_window}"
            "process.stdout.write(JSON.stringify(nearestWindow(105.5)));"
        )

        completed = subprocess.run(
            [node, "-e", runtime_source],
            input=json.dumps(embedded_data(html)["windows"]),
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(completed.stdout)[:2], [81, 130])

    def test_large_report_initializes_in_v8_without_argument_expansion(self):
        node = shutil.which("node")
        if node is None:
            self.skipTest("Node.js is unavailable for report JavaScript verification")
        window_count = 150_000
        windows = tuple(
            window(start, 0.95, 0.85, 1.0, 0.2)
            for start in range(1, window_count + 1)
        )
        html = self.render(windows=windows)
        data = embedded_data(html)
        self.assertEqual(len(data["windows"]), window_count)
        self.assertEqual(data["windows"][-1][:2], [150_000, 150_049])

        source = report_javascript(html)
        initialization = source[: source.index("let viewStart=")]
        runtime_source = (
            'const payload=require("fs").readFileSync(0,"utf8");'
            "const canvasStub={getContext:()=>({})};"
            "const document={getElementById:id=>id===\"geison-report-data\""
            "?{textContent:payload}:id===\"conservation-canvas\"?canvasStub:{}};"
            f"{initialization}"
            "process.stdout.write(String(genomeEnd));"
        )

        completed = subprocess.run(
            [node, "-e", runtime_source],
            input=embedded_data_source(html),
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(completed.stdout, "150049")

    def test_empty_report_has_an_accessible_visible_state(self):
        html = self.render(windows=())

        self.assertIn('id="empty-state"', html)
        self.assertIn('role="status"', html)
        self.assertIn('aria-live="polite"', html)
        self.assertNotIn('id="empty-state" hidden', html)
        self.assertEqual(embedded_data(html)["topWindows"], [])

    def test_browser_fixture_generator_runs_from_a_source_checkout(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_directory = Path(temporary_directory) / "fixture"
            completed = subprocess.run(
                [
                    sys.executable,
                    "integration_tests/generate_conservation_report.py",
                    str(output_directory),
                ],
                cwd=Path(__file__).resolve().parents[1],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            report_path = Path(completed.stdout.strip())
            self.assertEqual(report_path, (output_directory / "report.html").resolve())
            self.assertTrue(report_path.is_file())


if __name__ == "__main__":
    unittest.main()
