import pytest

from qpcr_pipeline.execution import (
    STAGE_DEPENDENCIES,
    STAGE_ORDER,
    ExecutionPolicy,
    plan_from_validity,
    required_reuse_boundary,
    transitive_descendants,
)


def _actions(decisions):
    return {decision.stage: decision.action for decision in decisions}


def _all_valid():
    return {stage: True for stage in STAGE_ORDER}


def test_panel_is_first_stage_and_input_depends_on_it():
    assert STAGE_ORDER[0] == "panel"
    assert STAGE_DEPENDENCIES["panel"] == ()
    assert STAGE_DEPENDENCIES["input"] == ("panel",)


def test_resume_invalid_panel_forces_every_downstream_stage():
    reusable = _all_valid()
    reusable["panel"] = False
    actions = _actions(plan_from_validity(ExecutionPolicy(resume=True), reusable))
    assert actions["panel"] == "RUN"
    assert all(actions[stage] == "FORCED" for stage in STAGE_ORDER[1:])


def test_specificity_descendants_only_include_ranking():
    assert transitive_descendants("specificity") == ("ranking",)


def test_inclusivity_force_does_not_force_specificity():
    assert set(transitive_descendants("inclusivity")) == {"ranking"}


def test_from_specificity_requires_inclusivity_branch():
    required = set(required_reuse_boundary("specificity"))
    assert "inclusivity" in required
    assert "primer_design" in required
    assert "qc" in required
    assert "specificity" not in required
    assert "ranking" not in required


def test_policy_rejects_resume_with_from_step():
    with pytest.raises(ValueError, match="resume.*from-step"):
        ExecutionPolicy(resume=True, from_step="alignment")


def test_policy_rejects_from_step_with_force_step():
    with pytest.raises(ValueError, match="from-step.*force-step"):
        ExecutionPolicy(from_step="alignment", force_step="alignment")


def test_policy_rejects_force_without_resume():
    with pytest.raises(ValueError, match="force-step.*resume"):
        ExecutionPolicy(force_step="alignment")


def test_policy_rejects_unknown_stage():
    with pytest.raises(ValueError, match="unknown pipeline stage"):
        ExecutionPolicy(resume=True, force_step="bogus")


def test_normal_run_runs_every_stage():
    decisions = plan_from_validity(ExecutionPolicy(), {})
    assert _actions(decisions) == {stage: "RUN" for stage in STAGE_ORDER}


def test_resume_reuses_all_valid_stages():
    decisions = plan_from_validity(ExecutionPolicy(resume=True), _all_valid())
    assert _actions(decisions) == {stage: "REUSE" for stage in STAGE_ORDER}


def test_resume_invalid_specificity_forces_only_specificity_and_ranking():
    reusable = _all_valid()
    reusable["specificity"] = False
    actions = _actions(plan_from_validity(ExecutionPolicy(resume=True), reusable))
    assert actions == {
        "panel": "REUSE",
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "REUSE",
        "alignment": "REUSE",
        "conservation": "REUSE",
        "primer_design": "REUSE",
        "inclusivity": "REUSE",
        "specificity": "RUN",
        "ranking": "FORCED",
    }


def test_resume_invalid_clustering_forces_entire_dependent_chain():
    reusable = _all_valid()
    reusable["clustering"] = False
    actions = _actions(plan_from_validity(ExecutionPolicy(resume=True), reusable))
    assert actions == {
        "panel": "REUSE",
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "RUN",
        "alignment": "FORCED",
        "conservation": "FORCED",
        "primer_design": "FORCED",
        "inclusivity": "FORCED",
        "specificity": "FORCED",
        "ranking": "FORCED",
    }


def test_resume_invalid_ranking_runs_only_ranking():
    reusable = _all_valid()
    reusable["ranking"] = False
    actions = _actions(plan_from_validity(ExecutionPolicy(resume=True), reusable))
    expected = {stage: "REUSE" for stage in STAGE_ORDER}
    expected["ranking"] = "RUN"
    assert actions == expected


def test_force_inclusivity_reuses_specificity_and_forces_ranking():
    policy = ExecutionPolicy(resume=True, force_step="inclusivity")
    actions = _actions(plan_from_validity(policy, _all_valid()))
    assert actions["inclusivity"] == "FORCED"
    assert actions["specificity"] == "REUSE"
    assert actions["ranking"] == "FORCED"
    assert actions["primer_design"] == "REUSE"


def test_from_step_specificity_runs_specificity_and_ranking_only():
    policy = ExecutionPolicy(from_step="specificity")
    actions = _actions(plan_from_validity(policy, _all_valid()))
    assert actions == {
        "panel": "REUSE",
        "input": "REUSE",
        "qc": "REUSE",
        "clustering": "REUSE",
        "alignment": "REUSE",
        "conservation": "REUSE",
        "primer_design": "REUSE",
        "inclusivity": "REUSE",
        "specificity": "FORCED",
        "ranking": "FORCED",
    }


def test_from_step_fails_before_execution_when_boundary_is_invalid():
    reusable = _all_valid()
    reusable["inclusivity"] = False
    with pytest.raises(ValueError, match="inclusivity"):
        plan_from_validity(ExecutionPolicy(from_step="specificity"), reusable)


def test_traversal_results_follow_stage_order():
    descendants = transitive_descendants("qc")
    positions = [STAGE_ORDER.index(stage) for stage in descendants]
    assert positions == sorted(positions)
