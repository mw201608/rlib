recoverGeneSymbolFromExcelDate=function(genes,species='human'){
#Recover Human gene symbols from excel date format
	Conv=c(`1-Dec`='DEC1', `2-Dec`='DEC2', `1-Mar`='MARCH1', `1-Sep`='SEPT1', `2-Mar`='MARCH2', `2-Sep`='SEPT2', `3-Mar`='MARCH3', `3-Sep`='SEPT3', `4-Mar`='MARCH4', `4-Sep`='SEPT4',
	`5-Mar`='MARCH5', `5-Sep`='SEPT5', `6-Mar`='MARCH6', `6-Sep`='SEPT6', `7-Mar`='MARCH7', `7-Sep`='SEPT7', `8-Mar`='MARCH8', `8-Sep`='SEPT8', `9-Mar`='MARCH9', `9-Sep`='SEPT9',
	`10-Mar`='MARCH10', `10-Sep`='SEPT10', `11-Mar`='MARCH11', `11-Sep`='SEPT11', `12-Sep`='SEPT12', `13-Sep`='SEPT13', `14-Sep`='SEPT14', `15-Sep`='SEPT15')
	if(species!='human') Conv=tools:::toTitleCase(tolower(Conv))
	id=genes %in% names(Conv)
	if(any(id)) genes[id]=Conv[genes[id]]
	genes
}
