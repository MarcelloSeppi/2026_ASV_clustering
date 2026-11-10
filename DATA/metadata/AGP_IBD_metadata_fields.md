# American Gut Project: Key Metadata Fields Related to IBD and Potential Confounders

This document lists ~50 metadata fields from the American Gut Project (AGP) that are most relevant for analyzing inflammatory bowel disease (IBD) and related confounding factors.  
The fields are grouped by thematic relevance.

---

## 🧫 Disease & Health Status (directly related or comorbid)

1. `ibd`  
2. `ibd_diagnosis`  
3. `ibs`  
4. `autoimmune`  
5. `acid_reflux`  
6. `diabetes`  
7. `thyroid`  
8. `cancer`  
9. `mental_illness`  
10. `liver_disease`  
11. `kidney_disease`  
12. `lung_disease`  
13. `cardiovascular_disease`  
14. `pm_health` (general self-reported health)

---

## 🍽️ Diet & Nutrition

15. `diet_type`  
16. `gluten`  
17. `lactose`  
18. `olive_oil`  
19. `whole_eggs`  
20. `red_meat_frequency`  
21. `vegetable_frequency`  
22. `fruit_frequency`  
23. `milk_cheese_frequency`  
24. `probiotic_frequency`  
25. `fermented_frequency` (or `fermented_consumed`)  
26. `alcohol_consumption`  
27. `vioscreen_calories`  
28. `vioscreen_protein`  
29. `vioscreen_fiber`  
30. `vioscreen_fat`

---

## ⚖️ Anthropometrics & Lifestyle

31. `bmi`  
32. `weight_kg`  
33. `height_cm`  
34. `age_years`  
35. `sex`  
36. `exercise_frequency`  
37. `sleep_quality`  
38. `sleep_duration`  
39. `smoking_frequency`  
40. `drinks_per_session`

---

## 🌍 Environmental & Socioeconomic Factors

41. `country`  
42. `state`  
43. `latitude`  
44. `longitude`  
45. `elevation`  
46. `env_biome`  
47. `env_feature`  
48. `setting` (urban/rural/etc.)  
49. `level_of_education`  
50. `economic_region`

---

## 💡 Notes on Inclusion Rationale

- **Direct disease relevance:** (1–14) capture IBD itself and conditions that often overlap immunologically or symptomatically.  
- **Dietary factors:** (15–30) reflect macro- and micro-nutritional inputs that are consistently linked to gut microbiome composition and inflammation.  
- **Body metrics & lifestyle:** (31–40) are major confounders for both microbiota and disease prevalence.  
- **Geographic & socioeconomic context:** (41–50) affect diet, exposures, and microbial ecology.

---

### Suggested Next Steps

- Rank fields by **predictive value or missingness** in the AGP dataset.  
- Use this list to design **confounder-controlled statistical models** for IBD-related microbiome analyses.
