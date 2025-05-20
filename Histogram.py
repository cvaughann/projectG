import matplotlib.pyplot as plt

# Your protein counts output dictionary.
protein_counts = {'HNRNPC': 242082, 'PTBP1': 116161, 'CELF2': 94643, 'AGO2': 80868, 'U2AF2': 72609, 'ELAVL1': 61774, 'CSTF2': 58553, 'UPF1': 57442, 'HNRNPA1': 55389, 'TARDBP': 29484, 'U2AF65': 25606, 'DDX3X': 25351, 'SRSF1': 24242, 'HNRNPH': 21802, 'CSTF2T': 17889, 'IGF2BP2': 14747, 'RBP_occupancy': 13865, 'FIP1L1': 13600, 'FUS': 12329, 'YTHDC1': 10635, 'LIN28B': 10165, 'NUDT21': 9763, 'RBM15B': 9723, 'CPSF7': 8215, 'ATXN2': 7918, 'CPSF6': 7887, 'YTHDF1': 5339, 'WDR33': 5164, 'TIAL1': 4719, 'ZC3H7B': 4539, 'MOV10': 4490, 'YTHDF2': 4434, 'DDX3': 4271, 'RBM15': 4156, 'ACI': 4068, 'DGCR8': 3760, 'QKI': 3568, 'ALYREF': 3223, 'HNRNPM': 3194, 'EWSR1': 3099, 'SRRM4': 3004, 'KHSRP': 2990, 'HNRNPD': 2960, 'CPSF3': 2947, 'CPSF1': 2906, 'TIA1': 2869, 'CHTOP': 2613, 'FMR1': 2556, 'HNRNPL': 2487, 'UTP3': 2393, 'MATR3': 2381, 'RBFOX1': 2377, 'TDP43': 2363, 'SAFB2': 2317, 'STAU1': 2316, 'SF3B1': 2236, 'DIS3L2': 2150, 'CPSF4': 2121, 'SUGP2': 2084, 'TAF15': 2023, 'RBM47': 2017, 'IGF2BP3': 1974, 'IGF2BP1': 1887, 'EIF4A3': 1807, 'HNRNPU': 1776, 'SLTM': 1471, 'AGGF1': 1436, 'CAPRIN1': 1436, 'EXOSC5': 1386, 'DROSHA': 1270, 'HLTF': 1117, 'TBRG4': 1082, 'RBFOX2': 1026, 'LIN28A': 983, 'LARP4': 968, 'GTF2F1': 922, 'NONO': 902, 'EIF3G': 896, 'MBNL2': 889, 'UCHL5': 881, 'RTCB': 857, 'FAM120A': 848, 'HNRNPK': 843, 'ILF3': 838, 'YTHDC2': 836, 'KHDRBS1': 803, 'ALKBH1': 777, 'PCBP2': 757, 'NIPBL': 738, 'DKC1': 703, 'XRN2': 690, 'PUM2': 686, 'SRSF7': 684, 'EIF3D': 674, 'PRPF4': 668, 'EFTUD2': 651, 'PPIL4': 648, 'SSB': 636, 'PRPF8': 617, 'BCLAF1': 615, 'EIF4G2': 613, 'AQR': 605, 'EIF3B': 592, 'SDOS': 571, 'LSM11': 566, 'SRSF10': 560, 'HNRNPUL1': 558, 'AARS': 549, 'SUPV3L1': 536, 'SRSF3': 534, 'DICER1': 527, 'GRWD1': 520, 'SFPQ': 520, 'CNBP': 499, 'TARBP2': 497, 'CPSF2': 485, 'PPIG': 480, 'NCBP3': 479, 'EIF3A': 475, 'NOLC1': 464, 'XPO5': 461, 'GEMIN5': 456, 'SAFB': 447, 'ALKBH5': 431, 'NCBP2': 424, 'NXF1': 422, 'NKRF': 413, 'FUBP3': 413, 'FASTKD2': 401, 'TRA2A': 392, 'METTL3': 385, 'SF3B4': 374, 'ZC3H11A': 373, 'HNRNPF': 373, 'YTHDF3': 372, 'HNRNPH1': 362, 'BUD13': 362, 'YBX3': 357, 'AKAP8L': 356, 'NOP58': 353, 'ZNF622': 352, 'PABPC4': 347, 'RBM22': 341, 'MKRN1': 334, 'FXR2': 334, 'DDX6': 333, 'DHX30': 332, 'FBL': 327, 'SRSF9': 326, 'WDR43': 323, 'ZRANB2': 320, 'U2AF1': 311, 'RBM10': 307, 'SMNDC1': 301, 'DDX24': 290, 'WRN': 288, 'NOP56': 267, 'DDX59': 254, 'METTL14': 253, 'DDX52': 243, 'TNRC6C': 232, 'PHF6': 227, 'XRCC6': 223, 'WDR3': 221, 'WTAP': 207, 'SMN': 203, 'TNRC6A': 202, 'YWHAG': 192, 'SF3A3': 192, 'AATF': 191, 'TNRC6B': 184, 'PABPN1': 181, 'RBPMS': 180, 'DDX51': 179, 'BCCIP': 177, 'UNR': 175, 'CPEB4': 167, 'DDX55': 161, 'FKBP4': 160, 'PCBP1': 159, 'SDAD1': 158, 'RBM5': 153, 'CDC40': 141, 'PRKRA': 140, 'FXR1': 137, 'TROVE2': 128, 'HNRNPA2B1': 125, 'EIF3H': 120, 'AKAP1': 117, 'PUM1': 111, 'ZNF800': 111, 'GPKOW': 106, 'RPS3': 101, 'LARP7': 98, 'SUB1': 92, 'FTO': 90, 'METAP2': 88, 'UTP18': 85, 'MTPAP': 68, 'NIP7': 64, 'NOL12': 56, 'STAU2': 53, 'NSUN2': 51, 'GRSF1': 50, 'APOBEC3C': 48, 'DDX21': 47, 'DDX42': 42, 'ZC3H8': 42, 'RPS11': 41, 'G3BP1': 34, 'SND1': 32, 'SERBP1': 32, 'GNL3': 31, 'IMP3': 29, 'PUS1': 28, 'NPM1': 28, 'ABCF1': 25, 'SLBP': 19, 'EZH2': 13, 'POLR2G': 9, 'SBDS': 7}

# ---------- Option 1: Bar Chart for Top 20 Proteins ----------

# Sort the dictionary by count, in descending order
sorted_counts = sorted(protein_counts.items(), key=lambda x: x[1], reverse=True)
proteins, counts = zip(*sorted_counts)

# Select the top 20 proteins
top_n = 20
top_proteins = proteins[:top_n]
top_counts = counts[:top_n]

plt.figure(figsize=(12, 6))
plt.bar(top_proteins, top_counts)
plt.xlabel('Protein')
plt.ylabel('Count')
plt.title('Top 20 Protein Counts')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
