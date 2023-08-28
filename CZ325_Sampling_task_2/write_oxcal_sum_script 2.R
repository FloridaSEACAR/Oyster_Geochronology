# DEP AGREEMENT NO. CZ325 Deliverable 2b (Radiocarbon sample size for age and time-averaging determination)
# helper function to create an OxCal program from a data from data frame of fractions and standard deviations
# df would have name, val, sd
write_oxcal_sum_script = function(df, path = ".", filename = "test",curvename = "bombBahamasto10000calBP_Marine20", 
  curvepath = here::here("CZ325_Sampling_task_2/bombBahamasto10000calBP_Marine20.14c"))
{
  
  if(!dir.exists(path)) dir.create(path)
  
  relative_name = paste0(paste(path,filename,sep = "/"),".oxcal")
  
  chunk1 <- ' 
Options()
 {
  Resolution=0.5;
  '
  
  Curve = paste0("Curve=",'"',curvepath,'"',";\n};")
  
chunk2 = '
 Plot()
 {
 '
chunk3 = sprintf('Curve("%s","%s");', curvename, curvepath)

chunk4 = '
  sum("Sum")
  {
  '
 
chunk5 = ""
for( ii in 1:nrow(df)){
  chunk5 = c(chunk5, sprintf('R_F14C("%s",%f,%f);\n', df[ii,1],df[ii,2], df[ii,3]))
  
}


chunk6 = '
  };
 };
'
  
code = c(chunk1, Curve, chunk2, chunk3, chunk4, chunk5, chunk6)
write(code, file = relative_name)
 return(code) 
}


