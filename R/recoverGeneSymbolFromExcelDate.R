recoverGeneSymbolFromExcelDate=function(genes, species='human'){
#Recover Human gene symbols from excel date format. This is not the best solution since 1-Mar can be either MARC1 (new official symbol MTARC1) or MARCH1 (new official symbol MARCHF1), while 2-Mar can be either MARC2 (new official symbol MTARC2) or MARCH2 (new official symbol MARCHF2)
	if(c('1-Mar') %in% genes) warning('1-Mar will be converted to MARCHF1. However it has two possibilities: MARC1 (new official symbol MTARC1) or MARCH1 (new official symbol MARCHF1). Please resolve it manually.\n')
	if(c('2-Mar') %in% genes) warning('2-Mar will be converted to MARCHF2. However it has two possibilities: MARC2 (new official symbol MTARC1) or MARCH2 (new official symbol MARCHF1). Please resolve it manually.\n')
	Conv=c(`1-Dec`='DELEC1', `2-Dec`='BHLHE41', `1-Mar`='MARCHF1', `1-Sep`='SEPTIN1', `2-Mar`='MARCHF2', `2-Sep`='SEPTIN2', `3-Mar`='MARCHF3', `3-Sep`='SEPTIN3', `4-Mar`='MARCHF4', `4-Sep`='SEPTIN4',
	`5-Mar`='MARCHF5', `5-Sep`='SEPTIN5', `6-Mar`='MARCHF6', `6-Sep`='SEPTIN6', `7-Mar`='MARCHF7', `7-Sep`='SEPTIN7', `8-Mar`='MARCHF8', `8-Sep`='SEPTIN8', `9-Mar`='MARCHF9', `9-Sep`='SEPTIN9',
	`10-Mar`='MARCHF10', `10-Sep`='SEPTIN10', `11-Mar`='MARCHF11', `11-Sep`='SEPTIN11', `12-Sep`='SEPTIN12', `13-Sep`='SEPTIN7P2', `14-Sep`='SEPTIN14', `15-Sep`='SELENOF')
	if(species!='human') Conv=to_title_case(tolower(Conv))
	id=genes %in% names(Conv)
	if(any(id)) genes[id]=Conv[genes[id]]
	genes
}
