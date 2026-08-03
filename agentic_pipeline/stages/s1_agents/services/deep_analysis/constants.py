"""Shared constants for the deep analysis pipeline."""

MAX_DISEASES_FOR_QUERY = 15

PENETRANCE_KEYWORDS = [
    "penetrance",
    '"complex trait"',
    '"complex disease"',
    "polygenic",
    "mendelian",
]

INHERITANCE_KEYWORDS = [
    "dominant",
    '"dominant negative"',
    "recessive",
    '"co-dominant"',
    "recessivity",
    "haploinsufficient",
    "haploinsufficiency",
]

MECHANISM_KEYWORDS = [
    '"gain-of-function"',
    '"loss-of-function"',
    '"dominant negative"',
    "gain",
    "loss",
]

ASSOCIATION_LABELS = {
    1: "Definitive evidence",
    2: "Sufficient evidence",
    3: "Moderate evidence",
    4: "Weak/conflicting evidence",
    5: "No disease found / default",
}

PENETRANCE_LEVEL_LABELS = {
    1: "Fully Mendelian penetrance",
    2: "High penetrance",
    3: "Moderate penetrance",
    4: "Low penetrance",
    5: "Complex trait / incomplete penetrance",
}

INHERITANCE_SCORE_LABELS = {
    1: "dominance",
    2: "incomplete dominance mostly dominant",
    3: "incomplete dominance",
    4: "co-dominant",
    5: "incomplete dominance mostly recessive",
    6: "recessive",
    7: "unknown-inheritance",
    8: "conflicting evidence",
}

ONSET_SCORE_LABELS = {
    1: "Prenatal",
    2: "Neonatal",
    3: "Infancy",
    4: "Childhood",
    5: "Adolescence",
    6: "Adulthood",
    7: "Late-onset",
    8: "Variable onset",
    9: "Unknown onset",
}

SEVERITY_LEVEL_LABELS = {
    1: "lethal or profoundly disabling clinical course",
    2: "severe but survivable presentation",
    3: "moderate morbidity",
    4: "mild manifestations",
    5: "very mild or subclinical impact",
}

INHERITANCE_PENALTIES = {
    1: 0,
    2: 1,
    3: 2,
    4: 2,
    5: 3,
    6: 4,
    7: 4,
    8: 4,
}

ONSET_PENALTIES = {
    1: 0,
    2: 0,
    3: 0,
    4: 1,
    5: 2,
    6: 5,
    7: 10,
    8: 5,
    9: 10,
}

ASSOCIATION_PENALTIES = {
    1: 0,
    2: 2,
    3: 5,
    4: 10,
}

DEFAULT_ALGO_LEVEL = 7
