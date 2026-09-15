# -*- coding: utf-8 -*-
"""Resolve {{KEY}} citations in MANUSCRIPT_NATMEMED.md into numbered refs.
- Only includes refs actually cited (no orphan padding).
- Prints compliance report: main-text words (NM definition), abstract words, ref count, display items.
"""
import io, re, sys, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
BASE = r"D:\南开CDC\task\1、监测殖民主义"

# ---------------- VERIFIED REFERENCE DICTIONARY (Nature style) ----------------
REFS = {
 "MEDEMA":   "Medema, G., Heijnen, L., Elsinga, G., Italiaander, R. & Brouwer, A. Presence of SARS-coronavirus-2 RNA in sewage and correlation with reported COVID-19 prevalence in the early stage of the epidemic. Environ. Sci. Technol. Lett. 7, 511\u2013516 (2020). https://doi.org/10.1021/acs.estlett.0c00357",
 "PECCIA":   "Peccia, J. et al. Measurement of SARS-CoV-2 RNA in wastewater tracks community infection dynamics. Nat. Biotechnol. 38, 1164\u20131167 (2020). https://doi.org/10.1038/s41587-020-0684-z",
 "KIRBY":    "Kirby, A. E. et al. Using wastewater surveillance data to support the COVID-19 response \u2014 United States, 2020\u20132021. MMWR Morb. Mortal. Wkly Rep. 70, 1242\u20131244 (2021). https://doi.org/10.15585/mmwr.mm7036a2",
 "HAN":      "Han, J. et al. Chinese urban wastewater surveillance system for early warning of infectious diseases: implementation and efficacy \u2014 January 2023\u2013June 2025. China CDC Wkly. 7, 1571\u20131576 (2025). https://doi.org/10.46234/ccdcw2025.267",
 "KESHAVIAH":"Keshaviah, A. et al. Wastewater monitoring can anchor global disease surveillance systems. Lancet Glob. Health 11, e976\u2013e981 (2023). https://doi.org/10.1016/S2214-109X(23)00170-5",
 "KNYAZEV":  "Knyazev, S. et al. Unlocking capacities of genomics for the COVID-19 response and future pandemics. Nat. Methods 19, 374\u2013380 (2022). https://doi.org/10.1038/s41592-022-01444-z",
 "NACHEGA":  "Nachega, J. B. et al. Advancing detection and response capacities for emerging and re-emerging pathogens in Africa. Lancet Infect. Dis. 23, e185\u2013e189 (2023). https://doi.org/10.1016/S1473-3099(22)00723-X",
 "MALLAPATY":"Mallapaty, S. Omicron-variant border bans ignore the evidence, say scientists. Nature 600, 199 (2021). https://doi.org/10.1038/d41586-021-03608-x",
 "DEOLIVEIRA":"de Oliveira, T. et al. Strengthening microbial genomics capacity in Africa for epidemic preparedness: key lessons from the COVID-19 pandemic. J. Glob. Health 16, 03012 (2026). https://doi.org/10.7189/jogh.16.03012",
 "PATINO":   "Pati\u00f1o, L. H. et al. Global and genetic diversity of SARS-CoV-2 in wastewater. Heliyon 10, e27452 (2024). https://doi.org/10.1016/j.heliyon.2024.e27452",
 "CARE":     "Carroll, S. R. et al. The CARE principles for Indigenous data governance. Data Sci. J. 19, 43 (2020). https://doi.org/10.5334/dsj-2020-043",
 "KROISS":   "Kroiss, S. J. et al. Assessing the sensitivity of the polio environmental surveillance system. PLoS One 13, e0208336 (2018). https://doi.org/10.1371/journal.pone.0208336",
 "SUAREZ":   "Suarez, C. et al. Detecting SARS-CoV-2 cryptic lineages using publicly available whole genome wastewater sequencing data. PLoS Pathog. 21, e1012850 (2025). https://doi.org/10.1371/journal.ppat.1012850",
 "SHU":      "Shu, Y. & McCauley, J. GISAID: global initiative on sharing all influenza data \u2014 from vision to reality. Eurosurveillance 22, 30494 (2017). https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494",
 "WHO_GENOMIC":"World Health Organization. Global genomic surveillance strategy for pathogens with pandemic and epidemic potential, 2022\u20132032. WHO, Geneva (2022). Available at: https://www.who.int/publications/i/item/9789240046979",
 "WHO_DASH": "World Health Organization. WHO COVID-19 dashboard: global data (cumulative deaths, as of 23 August 2026). https://data.who.int/dashboards/covid19/ (2026).",
 "UNWPP":    "United Nations Department of Economic and Social Affairs, Population Division. World Population Prospects 2022: Summary of Results. UN DESA/POP/2022/TR/NO.3, New York (2022). https://doi.org/10.18356/9789210014380",
 "ALKHOURI": "Alkhouri, N. B., Fisher, V. & Abuelezam, N. N. Equity in wastewater surveillance: assessment of the National Wastewater Surveillance System for COVID-19. AJPM Focus 5, 100434 (2026). https://doi.org/10.1016/j.focus.2025.100434",
 "RECKLING": "Reckling, S. et al. A geospatial analysis comparing wastewater-monitored sewershed and statewide populations for 32 states. PLoS Glob. Public Health 6, e0006243 (2026). https://doi.org/10.1371/journal.pgph.0006243",
 "AFRICACDC":"Makoni, M. Africa's $100-million pathogen genomics initiative. Lancet Microbe 1, e318 (2020). https://doi.org/10.1016/S2666-5247(20)30206-8",
 "JONES":    "Piasecki, J., Cheah, P. Y. & Arawi, T. Ownership of individual-level health data, data sharing, and data governance. BMC Med. Ethics 23, 104 (2022). https://doi.org/10.1186/s12910-022-00848-y",
 "EU2021":   "European Commission. Commission Implementing Decision (EU) 2021/953 of 14 June 2021. Off. J. Eur. Union L 212, 63\u201371 (2021).",
}

t = open(BASE + r"\MANUSCRIPT_NATMEMED.md", encoding="utf-8").read()
SENT = "\x00REFLIST\x00"
assert "{{REFLIST}}" in t
t = t.replace("## References\n\n{{REFLIST}}", "## References\n\n" + SENT)

# ---- numbering by first appearance ----
order, seen = [], set()
def repl(m):
    keys = re.findall(r"\{\{([A-Z0-9_]+)\}\}", m.group(0))
    nums = []
    for k in keys:
        if k not in REFS:
            raise RuntimeError("UNKNOWN KEY: " + k)
        if k not in seen:
            seen.add(k); order.append(k)
        nums.append(order.index(k) + 1)
    nums = sorted(set(nums))
    # compress consecutive runs
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j+1] == nums[j] + 1:
            j += 1
        out.append(str(nums[i]) if j == i else (f"{nums[i]},{nums[j]}" if j == i + 1 else f"{nums[i]}\u2013{nums[j]}"))
        i = j + 1
    return "^{" + ",".join(out) + "}"

body = re.sub(r"(?:\{\{[A-Z0-9_]+\}\})+", repl, t)
body = body.replace(SENT, "\n".join(f"{i+1}. {REFS[k]}" for i, k in enumerate(order)))

open(BASE + r"\NM_submission_final.md", "w", encoding="utf-8").write(body)

# ---- compliance report ----
def wc(s):
    s = re.sub(r"\^\{[^}]*\}", " ", s)
    s = re.sub(r"\{\{[^}]+\}\}", " ", s)
    return len(re.findall(r"[A-Za-z][A-Za-z'\-]*", s))

main = body.split("## Introduction")[1].split("## Online Methods")[0]
main = "Introduction" + main
abstract = body.split("## Abstract")[1].split("## Introduction")[0]
refs = body.split("## References")[1].split("## Figure legends")[0]
ref_lines = [l for l in refs.splitlines() if re.match(r"^\d+\. ", l)]

print("=" * 62)
print(f"Main text words (intro+results+discussion): {wc(main):,}   [NM limit 4,000]")
print(f"Abstract words:                             {wc(abstract)}   [NM limit 150]")
print(f"References:                                 {len(ref_lines)}   [NM guideline ~60]")
print(f"Placeholders left: {body.count('{{')}   [must be 0]")
print("=" * 62)
cited = set(order)
orphan = [k for k in REFS if k not in cited]
print("Refs defined but never cited:", orphan)
nums = [int(re.match(r"^(\d+)\. ", l).group(1)) for l in ref_lines]
print("Ref numbering contiguous:", nums == list(range(1, len(nums) + 1)))
print("Citation first-appearance order:", order)
