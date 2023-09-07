# https://github.com/broadinstitute/tensorqtl


#module load lang/Python/3.9.5-GCCcore-10.3.0
#module load lang/R/4.1.0-foss-2021a

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

# define paths to data
plink_prefix_path = '/PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf5_EUR'
expression_bed = '/PROJECTES/BISC_OMICS/analyses/BiSC_21/006_PlamiQTLs_AG/data/BISC_miRNA.bed.gz'
covariates_file = '/PROJECTES/BISC_OMICS/analyses/BiSC_21/006_PlamiQTLs_AG/data/BISC_Cells_covariates.txt' #"b_sexo", "gestage_0y_c_days", "Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC", "PC1", "PC2", "PC3", "PC4", "PC5"
prefix = 'BISC_EUR_mat_miQTL'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

variant_df['chrom'] = variant_df['chrom'].astype('O')
 
list_sample = ['BISC_0682','BISC_0687','BISC_0700','BISC_0708','BISC_0713','BISC_0723','BISC_0725','BISC_0727','BISC_0728','BISC_0734','BISC_0745','BISC_0759','BISC_0765','BISC_0767','BISC_0769','BISC_0776','BISC_0781','BISC_0782','BISC_0792','BISC_0796','BISC_0798','BISC_0802','BISC_0805','BISC_0807','BISC_0811','BISC_0815','BISC_0820','BISC_0821','BISC_0822','BISC_0826','BISC_0827_bis','BISC_0833','BISC_0842_bis','BISC_0848_bis','BISC_0849','BISC_0853','BISC_0856','BISC_0859','BISC_0862','BISC_0868','BISC_0869','BISC_0870','BISC_0875','BISC_0877','BISC_0879','BISC_0884','BISC_0893','BISC_0898_bis','BISC_0904','BISC_0907','BISC_0908','BISC_0909','BISC_0910','BISC_0911_bis','BISC_0916','BISC_0922','BISC_0923','BISC_0935','BISC_0937','BISC_0938','BISC_0948','BISC_0951_bis','BISC_0954','BISC_0956','BISC_0957','BISC_0961','BISC_0962','BISC_0963_bis','BISC_0965','BISC_0970','BISC_0971','BISC_0972','BISC_0977','BISC_0983','BISC_0985','BISC_0988','BISC_1001','BISC_1003','BISC_1004','BISC_1006','BISC_1007','BISC_1009','BISC_1010_bis','BISC_1013_bis','BISC_1015','BISC_1017','BISC_1018','BISC_1027','BISC_1032','BISC_1033','BISC_1036','BISC_1037','BISC_1043','BISC_1044','BISC_1046','BISC_1050','BISC_1056','BISC_1058','BISC_1065','BISC_1066','BISC_1067','BISC_1072','BISC_1076','BISC_1082','BISC_1085','BISC_1087','BISC_1088','BISC_1089','BISC_1094','BISC_1095','BISC_1103','BISC_1106','BISC_1109','BISC_1110','BISC_1111','BISC_1115','BISC_1117','BISC_1118','BISC_1122','BISC_1127','BISC_1131','BISC_1138','BISC_1139','BISC_1153','BISC_1159','BISC_1164','BISC_1168','BISC_1169','BISC_1171','BISC_1182','BISC_1186','BISC_1190','BISC_1193','BISC_1201','BISC_1202','BISC_1204','BISC_1209','BISC_1214','BISC_1217','BISC_1224','BISC_1226','BISC_1229','BISC_1231_bis','BISC_1236','BISC_1237','BISC_1238_bis','BISC_1243','BISC_1246','BISC_1247','BISC_1255','BISC_1268','BISC_1275','BISC_1279','BISC_1284','BISC_1287','BISC_1288_bis','BISC_1289_bis','BISC_1292','BISC_1294','BISC_1297','BISC_1298_bis','BISC_1300','BISC_1308_bis','BISC_1309','BISC_1312','BISC_1315_bis','BISC_1327','BISC_1330','BISC_1332','BISC_1336','BISC_1340','BISC_1350','BISC_1353','BISC_1354','BISC_1357','BISC_1367','BISC_1368','BISC_1371','BISC_1379','BISC_1381','BISC_1385','BISC_1386','BISC_1392','BISC_1393','BISC_1396','BISC_1397','BISC_1398','BISC_1410','BISC_1413','BISC_1416','BISC_1418','BISC_1419','BISC_1429','BISC_1440','BISC_1447','BISC_1454','BISC_1457','BISC_1459','BISC_1460','BISC_1464','BISC_1473','BISC_1476','BISC_1477','BISC_1482','BISC_1486','BISC_1488_bis','BISC_1494','BISC_1500','BISC_1508','BISC_1509_bis','BISC_1510','BISC_1513','BISC_1522','BISC_1523','BISC_1532','BISC_1533','BISC_1534','BISC_1537','BISC_1538','BISC_1545','BISC_1550_bis','BISC_1551','BISC_1553','BISC_1554','BISC_1561','BISC_1563','BISC_1574','BISC_1575','BISC_1577','BISC_1581','BISC_1582','BISC_1585','BISC_1592','BISC_1596','BISC_1597','BISC_1603','BISC_1608','BISC_1611','BISC_1614','BISC_1621','BISC_1622','BISC_1623','BISC_1624','BISC_1629','BISC_1631','BISC_1636','BISC_1644','BISC_1652','BISC_1653','BISC_1655','BISC_1656','BISC_1657','BISC_1658','BISC_1661','BISC_1662','BISC_1666','BISC_1673','BISC_1676','BISC_1677','BISC_1682','BISC_1684','BISC_1689','BISC_1690','BISC_1691','BISC_1692','BISC_1694','BISC_1695','BISC_1700','BISC_1702','BISC_1707','BISC_1710']

 
genotype_df2 = genotype_df[list_sample]
covariates_df2 = covariates_df.loc[list_sample,:]
phenotype_df2 = phenotype_df[list_sample]

#trans-QTL mapping
trans_df = trans.map_trans(genotype_df2, phenotype_df2, covariates_df2, batch_size=10000,
    return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)
    
trans_df.to_csv('bisc_m_TensorQTL_pvalue1E-05.csv', index=False)






