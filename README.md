# Wet-Lab Agent Skills

Agent Skills that take routine bench results from instrument output to checked statistics, a figure, and a
draft report. They run offline, and each one is built to prevent the mistakes its analysis commonly attracts.

| Skill | Takes in | Gives you |
| --- | --- | --- |
| [qpcr-analysis](skills/qpcr-analysis/SKILL.md) | QuantStudio, CFX, or similar exports | Efficiency, reference ranking, fold change with intervals, plate QC |
| [elisa-analysis](skills/elisa-analysis/SKILL.md) | Plate reader grid and layout | 4PL or 5PL curve, working range, concentrations with intervals, censoring |
| [western-blot-densitometry](skills/western-blot-densitometry/SKILL.md) | Band density tables | Within-membrane normalization, linearity, saturation and loading QC |
| [candidate-gene-association](skills/candidate-gene-association/SKILL.md) | Genotype tables or counts | Power, Hardy-Weinberg, association models, meta-analysis, forest plots |

## Install

Copy any folder from `skills/` into your agent's skills directory. For Claude Code that is `~/.claude/skills/`.
Then ask in plain language, for example "analyse my qPCR export".

Requirements are Python 3.11 or newer, plus `scipy` for qPCR, Western blot, and ELISA, and `numpy` for ELISA.
`matplotlib` is needed for figures and `openpyxl` for .xlsx import.

## Status

The statistics are tested against independent values, and the tests were themselves checked by injecting
bugs. **The importers follow documented export layouts and synthetic examples, and have not yet been tried on
real instrument files.** Check your first run against the original data, and treat every report as a draft
for a person to review. Thresholds are working defaults unless a reference says otherwise. These skills are
research tools and are not validated for clinical or regulated use.

## Develop

```bash
python -m pip install pytest scipy numpy matplotlib openpyxl
python tests/run_all.py
```

Structure follows the [Agent Skills specification](https://agentskills.io/specification), and conventions
follow [K-Dense scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills).

---

# مهارات المختبر للوكلاء

مهارات تحوّل نتائج المختبر اليومية من ملف الجهاز إلى إحصاء مفحوص ورسم بياني ومسودة تقرير، وتشتغل بدون إنترنت.

انسخ أي مجلد من `skills/` إلى مجلد مهارات الوكيل، وفي Claude Code هو `~/.claude/skills/`، وبعدها اطلب بكلامك العادي،
مثل: «حلل لي ملف qPCR». اختبارات الاستيراد تمت على أمثلة مصطنعة، فراجع أول نتيجة مع ملفك الأصلي، واعتبر كل تقرير مسودة.

## License

MIT. See [LICENSE](LICENSE).
