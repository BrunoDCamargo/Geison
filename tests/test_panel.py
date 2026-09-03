from dataclasses import replace

import pytest

from qpcr_pipeline.panel import (
    DiagnosticContext,
    PanelDefinition,
    PanelNonTarget,
    PanelTarget,
    SequenceSelectionProvenance,
    TargetGroup,
    validate_panel_definition,
)


def valid_panel() -> PanelDefinition:
    selection = SequenceSelectionProvenance(
        dataset_role="DESIGN",
        method="manual_fixture",
        source="unit-test",
        details=("representative seed",),
    )
    return PanelDefinition(
        target=PanelTarget(
            name="Target virus",
            taxid=1001,
            mode="broad_detection",
            subtype=None,
            groups=(
                TargetGroup(
                    name="group-a",
                    required=True,
                    dataset_roles=("DESIGN", "CHALLENGE"),
                    reasons=("target_diversity",),
                    proposed_by=("manual",),
                    sequence_selection=(selection,),
                ),
            ),
        ),
        non_targets=(
            PanelNonTarget(
                name="Neighbor virus",
                taxid=2001,
                criticality="CRITICAL",
                dataset_roles=("DESIGN", "CHALLENGE"),
                reasons=("phylogenetic_neighbor",),
                proposed_by=("manual",),
                sequence_selection=(),
            ),
        ),
        diagnostic_context=DiagnosticContext(
            syndrome="febrile illness",
            geography="test-region",
            sample_type="serum",
            vector="mosquito",
        ),
    )


def test_valid_panel_definition_is_accepted():
    validate_panel_definition(valid_panel())


def test_subtype_specific_requires_subtype():
    panel = valid_panel()
    invalid = replace(
        panel,
        target=replace(panel.target, mode="subtype_specific", subtype=None),
    )
    with pytest.raises(ValueError, match="subtype_specific.*subtype"):
        validate_panel_definition(invalid)


def test_broad_detection_rejects_subtype():
    panel = valid_panel()
    invalid = replace(panel, target=replace(panel.target, subtype="group-a"))
    with pytest.raises(ValueError, match="broad_detection.*subtype"):
        validate_panel_definition(invalid)


def test_non_target_cannot_duplicate_target_taxid():
    panel = valid_panel()
    duplicate = replace(panel.non_targets[0], name="Target alias", taxid=1001)
    invalid = replace(panel, non_targets=(duplicate,))
    with pytest.raises(ValueError, match="target.*non-target.*TaxID"):
        validate_panel_definition(invalid)


def test_dataset_roles_must_be_unique():
    panel = valid_panel()
    invalid_group = replace(
        panel.target.groups[0],
        dataset_roles=("DESIGN", "DESIGN"),
    )
    invalid = replace(
        panel,
        target=replace(panel.target, groups=(invalid_group,)),
    )
    with pytest.raises(ValueError, match="dataset_roles.*unique"):
        validate_panel_definition(invalid)


@pytest.mark.parametrize("name", ["", "   "])
def test_target_name_must_be_non_blank(name):
    panel = valid_panel()
    with pytest.raises(ValueError, match="target name.*non-blank"):
        validate_panel_definition(
            replace(panel, target=replace(panel.target, name=name))
        )


def test_target_group_name_must_be_non_blank():
    panel = valid_panel()
    group = replace(panel.target.groups[0], name=" ")
    with pytest.raises(ValueError, match="target group 1 name.*non-blank"):
        validate_panel_definition(
            replace(panel, target=replace(panel.target, groups=(group,)))
        )


def test_non_target_name_must_be_non_blank():
    panel = valid_panel()
    item = replace(panel.non_targets[0], name=" ")
    with pytest.raises(ValueError, match="non-target 1 name.*non-blank"):
        validate_panel_definition(replace(panel, non_targets=(item,)))


@pytest.mark.parametrize("taxid", [0, -1, True])
def test_target_taxid_must_be_a_positive_integer(taxid):
    panel = valid_panel()
    with pytest.raises(ValueError, match="target TaxID.*positive integer"):
        validate_panel_definition(
            replace(panel, target=replace(panel.target, taxid=taxid))
        )


@pytest.mark.parametrize("taxid", [0, -1, False])
def test_non_target_taxid_must_be_a_positive_integer(taxid):
    panel = valid_panel()
    item = replace(panel.non_targets[0], taxid=taxid)
    with pytest.raises(ValueError, match="non-target 1 TaxID.*positive integer"):
        validate_panel_definition(replace(panel, non_targets=(item,)))


def test_target_group_names_are_unique_after_casefolding():
    panel = valid_panel()
    duplicate = replace(panel.target.groups[0], name=" GROUP-A ")
    target = replace(
        panel.target,
        groups=(panel.target.groups[0], duplicate),
    )
    with pytest.raises(ValueError, match="target group names.*unique"):
        validate_panel_definition(replace(panel, target=target))


def test_non_target_names_are_unique_after_casefolding():
    panel = valid_panel()
    duplicate = replace(panel.non_targets[0], name=" NEIGHBOR VIRUS ", taxid=2002)
    with pytest.raises(ValueError, match="non-target names.*unique"):
        validate_panel_definition(
            replace(panel, non_targets=(panel.non_targets[0], duplicate))
        )


def test_non_target_taxids_must_be_unique_when_present():
    panel = valid_panel()
    duplicate = replace(panel.non_targets[0], name="Other virus")
    with pytest.raises(ValueError, match="non-target TaxIDs.*unique"):
        validate_panel_definition(
            replace(panel, non_targets=(panel.non_targets[0], duplicate))
        )


def test_non_target_cannot_duplicate_target_name_after_casefolding():
    panel = valid_panel()
    duplicate = replace(panel.non_targets[0], name=" TARGET VIRUS ")
    with pytest.raises(ValueError, match="target.*non-target names.*distinct"):
        validate_panel_definition(replace(panel, non_targets=(duplicate,)))


def test_invalid_target_mode_is_rejected():
    panel = valid_panel()
    with pytest.raises(ValueError, match="target mode.*unsupported"):
        validate_panel_definition(
            replace(panel, target=replace(panel.target, mode="other"))
        )


def test_invalid_criticality_is_rejected():
    panel = valid_panel()
    item = replace(panel.non_targets[0], criticality="OTHER")
    with pytest.raises(ValueError, match="criticality.*unsupported"):
        validate_panel_definition(replace(panel, non_targets=(item,)))


@pytest.mark.parametrize("owner", ["target_group", "non_target"])
def test_invalid_dataset_role_is_rejected(owner):
    panel = valid_panel()
    if owner == "target_group":
        group = replace(panel.target.groups[0], dataset_roles=("TRAIN",))
        invalid = replace(panel, target=replace(panel.target, groups=(group,)))
    else:
        item = replace(panel.non_targets[0], dataset_roles=("TRAIN",))
        invalid = replace(panel, non_targets=(item,))
    with pytest.raises(ValueError, match="dataset_roles.*unsupported"):
        validate_panel_definition(invalid)


@pytest.mark.parametrize("owner", ["target_group", "non_target"])
def test_dataset_roles_must_not_be_empty(owner):
    panel = valid_panel()
    if owner == "target_group":
        group = replace(panel.target.groups[0], dataset_roles=())
        invalid = replace(panel, target=replace(panel.target, groups=(group,)))
    else:
        item = replace(panel.non_targets[0], dataset_roles=())
        invalid = replace(panel, non_targets=(item,))
    with pytest.raises(ValueError, match="dataset_roles.*must not be empty"):
        validate_panel_definition(invalid)


@pytest.mark.parametrize("field", ["reasons", "proposed_by"])
@pytest.mark.parametrize("owner", ["target_group", "non_target"])
def test_reason_and_proposer_values_must_be_non_blank(field, owner):
    panel = valid_panel()
    if owner == "target_group":
        group = replace(panel.target.groups[0], **{field: (" ",)})
        invalid = replace(panel, target=replace(panel.target, groups=(group,)))
    else:
        item = replace(panel.non_targets[0], **{field: (" ",)})
        invalid = replace(panel, non_targets=(item,))
    with pytest.raises(ValueError, match=f"{field}.*non-blank"):
        validate_panel_definition(invalid)


@pytest.mark.parametrize("field", ["syndrome", "geography", "sample_type", "vector"])
def test_optional_diagnostic_context_values_must_be_non_blank(field):
    panel = valid_panel()
    context = replace(panel.diagnostic_context, **{field: " "})
    with pytest.raises(ValueError, match=f"diagnostic_context.{field}.*non-blank"):
        validate_panel_definition(replace(panel, diagnostic_context=context))


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"dataset_role": "TRAIN"}, "dataset_role.*unsupported"),
        ({"method": " "}, "method.*non-blank"),
        ({"source": " "}, "source.*non-blank"),
        ({"details": (" ",)}, "details.*non-blank"),
        ({"details": ["seed"]}, "details.*tuple"),
    ],
)
def test_sequence_selection_provenance_fields_are_validated(changes, message):
    panel = valid_panel()
    selection = replace(panel.target.groups[0].sequence_selection[0], **changes)
    group = replace(panel.target.groups[0], sequence_selection=(selection,))
    invalid = replace(panel, target=replace(panel.target, groups=(group,)))
    with pytest.raises(ValueError, match=message):
        validate_panel_definition(invalid)


def test_required_must_be_exactly_boolean():
    panel = valid_panel()
    group = replace(panel.target.groups[0], required=1)
    with pytest.raises(ValueError, match="required.*boolean"):
        validate_panel_definition(
            replace(panel, target=replace(panel.target, groups=(group,)))
        )


def test_target_groups_must_not_be_empty():
    panel = valid_panel()
    with pytest.raises(ValueError, match="target groups.*non-empty tuple"):
        validate_panel_definition(
            replace(panel, target=replace(panel.target, groups=()))
        )
