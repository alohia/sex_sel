bys caste: egen caste_tot = sum(t)
egen temp = sum(t)
gen caste_weight = caste_tot/temp
drop temp
gen wt_PG = d_girl6 * caste_weight 
collapse (sum) wt_PG , by(inc_cls)
