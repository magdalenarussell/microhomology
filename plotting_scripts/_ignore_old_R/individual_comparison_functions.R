map_scanID_to_localID <- function(scanIDs_to_convert){
    # This function converts subject scanIDs to localIDs
    ID_map_file = fread(ID_MAPPING_FILE)
    require(plyr)
    converted_IDs = plyr::mapvalues(scanIDs_to_convert, ID_map_file$scanID, ID_map_file$localID)
    return(converted_IDs)
}

compile_all_genotypes_snp_list <- function(snp_list){
    require(SNPRelate)
    require(gdsfmt)
    snp_gds_file = openfn.gds(SNP_GDS_FILE)
    genotypes = snpgdsGetGeno(snp_gds_file, snp.id = snp_list, with.id = TRUE)
    closefn.gds(snp_gds_file)

    genotype_matrix = genotypes$genotype
    rownames(genotype_matrix) = map_scanID_to_localID(genotypes$sample.id)
    colnames(genotype_matrix) = genotypes$snp.id
    genotype_dt = data.table(localID = row.names(genotype_matrix), genotype_matrix)
    return(genotype_dt)
}
