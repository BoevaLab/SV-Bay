import glob
#files=['sampe_1.sam','sampe_2_.sam','sampe_3_.sam','sampe_4_.sam','sampe_6_.sam','sampe_7_.sam','sampe_8_.sam','081224_EAS188_0073_FC30TTCAAXX_1_new.sam','081224_EAS188_0073_FC30TTCAAXX_3_new.sam']
chrom_names=[]
#for i in range(1,23,1):
#	chrom_names.append('chr'+str(i))

#chrom_names.append('chrY')
#chrom_names.append('chrX')
#chrom_names.remove('chr20')
#chrom_names.remove('chr7')
#chrom_names =['chr20']
#print chrom_names
chrom_names = ['NODE_101_length_14188_cov_4.233930.fa', 'NODE_179_length_9479_cov_3.446671.fa', 'NODE_250_length_7033_cov_4.151287.fa', 'NODE_104_length_5184_cov_3.857446.fa', 'NODE_17_length_1080_cov_3.637963.fa', 'NODE_252_length_2736_cov_3.790936.fa', 'NODE_106_length_3199_cov_4.113473.fa', 'NODE_182_length_6780_cov_3.785988.fa', 'NODE_255_length_4437_cov_3.717827.fa', 'NODE_107_length_3905_cov_4.000000.fa', 'NODE_187_length_15556_cov_3.934237.fa', 'NODE_256_length_13513_cov_3.902316.fa', 'NODE_108_length_3400_cov_3.532647.fa', 'NODE_188_length_4981_cov_3.468380.fa', 'NODE_258_length_5201_cov_3.787349.fa', 'NODE_110_length_5089_cov_4.254863.fa', 'NODE_18_length_2280_cov_4.506579.fa', 'NODE_25_length_6296_cov_4.088310.fa', 'NODE_113_length_7208_cov_3.724334.fa', 'NODE_190_length_16726_cov_3.944338.fa', 'NODE_261_length_3732_cov_3.874866.fa', 'NODE_114_length_6840_cov_4.206433.fa', 'NODE_192_length_3773_cov_3.499603.fa', 'NODE_263_length_2881_cov_3.549809.fa', 'NODE_115_length_2449_cov_3.543895.fa', 'NODE_193_length_1949_cov_3.398666.fa', 'NODE_265_length_6579_cov_3.935857.fa', 'NODE_116_length_5490_cov_3.986521.fa', 'NODE_194_length_21946_cov_3.779869.fa', 'NODE_268_length_6040_cov_3.714570.fa', 'NODE_117_length_7944_cov_4.199647.fa', 'NODE_195_length_8677_cov_3.738734.fa', 'NODE_272_length_1616_cov_3.925124.fa', 'NODE_120_length_1808_cov_3.623894.fa', 'NODE_197_length_11657_cov_3.711161.fa', 'NODE_276_length_2039_cov_4.203531.fa', 'NODE_123_length_3625_cov_4.095724.fa', 'NODE_198_length_12557_cov_4.068089.fa', 'NODE_279_length_1545_cov_3.806473.fa', 'NODE_124_length_15482_cov_4.010529.fa', 'NODE_200_length_8502_cov_3.978711.fa', 'NODE_283_length_4163_cov_3.801826.fa', 'NODE_125_length_1085_cov_3.274654.fa', 'NODE_201_length_22760_cov_3.767311.fa', 'NODE_286_length_1562_cov_3.814341.fa', 'NODE_127_length_2906_cov_3.872677.fa', 'NODE_202_length_2151_cov_4.605300.fa', 'NODE_292_length_10666_cov_3.957716.fa', 'NODE_128_length_2511_cov_3.477897.fa', 'NODE_204_length_8292_cov_3.959841.fa', 'NODE_293_length_30040_cov_3.986119.fa', 'NODE_130_length_7360_cov_3.967391.fa', 'NODE_206_length_19598_cov_3.938718.fa', 'NODE_31_length_1236_cov_3.881877.fa', 'NODE_133_length_3213_cov_3.886399.fa', 'NODE_210_length_13118_cov_4.030111.fa', 'NODE_33_length_1944_cov_4.741769.fa', 'NODE_135_length_7663_cov_3.897429.fa', 'NODE_211_length_13153_cov_3.510226.fa', 'NODE_35_length_1841_cov_4.451385.fa', 'NODE_137_length_3537_cov_3.956743.fa', 'NODE_212_length_1439_cov_3.748436.fa', 'NODE_36_length_1245_cov_5.691566.fa', 'NODE_138_length_6411_cov_3.954453.fa', 'NODE_213_length_11016_cov_3.904503.fa', 'NODE_38_length_2075_cov_3.667952.fa', 'NODE_13_length_1039_cov_3.559191.fa', 'NODE_215_length_13407_cov_3.911017.fa', 'NODE_48_length_2200_cov_3.869091.fa', 'NODE_140_length_2211_cov_4.131162.fa', 'NODE_217_length_8074_cov_3.765791.fa', 'NODE_49_length_1394_cov_4.038737.fa', 'NODE_144_length_2576_cov_4.541537.fa', 'NODE_219_length_15920_cov_3.584234.fa', 'NODE_50_length_1812_cov_3.320088.fa', 'NODE_145_length_5865_cov_4.159591.fa', 'NODE_21_length_1378_cov_3.547170.fa', 'NODE_51_length_1301_cov_4.241353.fa', 'NODE_146_length_5343_cov_3.642710.fa', 'NODE_220_length_2541_cov_3.456513.fa', 'NODE_52_length_2777_cov_3.598488.fa', 'NODE_147_length_2112_cov_4.052557.fa', 'NODE_223_length_2622_cov_3.490465.fa', 'NODE_53_length_3956_cov_3.967897.fa', 'NODE_148_length_4498_cov_3.676523.fa', 'NODE_225_length_11656_cov_4.006692.fa', 'NODE_56_length_6765_cov_3.903326.fa', 'NODE_150_length_21950_cov_3.830433.fa', 'NODE_226_length_2294_cov_4.039669.fa', 'NODE_57_length_1229_cov_5.803905.fa', 'NODE_155_length_3398_cov_3.577693.fa', 'NODE_228_length_16792_cov_3.992020.fa', 'NODE_66_length_1372_cov_3.696793.fa', 'NODE_157_length_3189_cov_3.599247.fa', 'NODE_229_length_6655_cov_3.675282.fa', 'NODE_67_length_1754_cov_3.173318.fa', 'NODE_159_length_10054_cov_4.023076.fa', 'NODE_231_length_19195_cov_3.853816.fa', 'NODE_72_length_1116_cov_2.767025.fa', 'NODE_160_length_9937_cov_3.930562.fa', 'NODE_232_length_7373_cov_3.636647.fa', 'NODE_73_length_3813_cov_3.444007.fa', 'NODE_161_length_10937_cov_3.858279.fa', 'NODE_233_length_6170_cov_3.877958.fa', 'NODE_7_length_3652_cov_4.288883.fa', 'NODE_163_length_18198_cov_3.827564.fa', 'NODE_234_length_24373_cov_3.784352.fa', 'NODE_80_length_5278_cov_4.184729.fa', 'NODE_165_length_11268_cov_3.766596.fa', 'NODE_235_length_6281_cov_4.256010.fa', 'NODE_87_length_9569_cov_3.779078.fa', 'NODE_166_length_4273_cov_3.088931.fa', 'NODE_237_length_2879_cov_3.267801.fa', 'NODE_8_length_1243_cov_5.734513.fa', 'NODE_167_length_7819_cov_4.150531.fa', 'NODE_238_length_25248_cov_3.834799.fa', 'NODE_91_length_10939_cov_3.958132.fa', 'NODE_168_length_1375_cov_3.506909.fa', 'NODE_242_length_15521_cov_4.142452.fa', 'NODE_93_length_1373_cov_3.665696.fa', 'NODE_170_length_14269_cov_4.055295.fa', 'NODE_245_length_14977_cov_4.291046.fa', 'NODE_94_length_4752_cov_4.045455.fa', 'NODE_176_length_4231_cov_4.103049.fa', 'NODE_246_length_10015_cov_3.979830.fa', 'NODE_96_length_1036_cov_4.551158.fa', 'NODE_178_length_10768_cov_3.818722.fa', 'NODE_248_length_32660_cov_4.043356.fa']
data_files = glob.glob('/Users/daria/Dropbox/data/sam_files/simulate/' + '*.sam')
#data_files = glob.glob('/home/dasha/Boeva/data/sam_files/files_sam_chrom/' + '*7.sam')
for chrom in chrom_names:
	print 'Now for this chromosome = ', chrom
	new_file = open('/Users/dasha/Dropbox/GCAss_DB/sample1/'+chrom+'.sam','w')
	for file_curr in data_files:
		print 'Now for this file = ', file_curr
		flag=0
		j=0
		n=1
		for line in open(file_curr,'r'):
			j+=1
			if line[0]=='@':
				#print 'I write!!!'
				new_file.write(line)
			elif not flag:
				line1=line
				flag=1
			elif flag:
				line2=line
				temp1=line1.split()
				temp2=line2.split()
				if temp2[2]== chrom or temp1[2]== chrom:
					new_file.write(line1)
					new_file.write(line2)
				flag=0
			if j == 100000*n:
				print 'prossed ', 100000*n,' lines'
				n+=1	
	new_file.close()			
print chrom_names
