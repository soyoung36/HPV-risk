# Processed calibration data

`proc_KoreaHPVcancer_1999-2022.csv` contains annual counts for 1999-2022.
Columns are calendar years. The rows are:

- `Year`: calendar year.
- `fIncidence`: annual cervical cancer cases.
- `fDeath`: annual cervical cancer deaths.
- `mIncidence`: processed annual HPV-attributable male cancer cases (noninteger).

Data sources:

1. Ministry of Health and Welfare, Cancer Registration Statistics:
   https://kosis.kr/statHtml/statHtml.do?orgId=117&tblId=DT_117N_A00025&conn_path=I2
   Accessed 16 January 2026.
2. Statistics Korea, Cause of Death Statistics:
   https://kosis.kr/statHtml/statHtml.do?orgId=101&tblId=DT_1B34E01&conn_path=I2
   Accessed 16 January 2026.

The male series is the annual sum of site-specific male incidence multiplied
by the following HPV-attributable fractions:

`0.88*anus(C21) + 0.50*penis(C60) + 0.308*oropharynx(C01,C09-C10)`
`+ 0.024*larynx(C32) + 0.022*oral_cavity(C02-C06)`.

Fractions: de Martel C, Plummer M, Vignat J, Franceschi S (2017),
*International Journal of Cancer* 141, 664-670.
https://doi.org/10.1002/ijc.30716
