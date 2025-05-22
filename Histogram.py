import matplotlib.pyplot as plt

# Your protein counts output dictionary.
protein_counts = {'AGO2': 3794, 'HNRNPC': 3627, 'UPF1': 2097, 'U2AF65': 1894, 'ELAVL1': 1736, 'U2AF2': 1603, 'SRSF1': 1277, 'CSTF2': 1163, 'HNRNPH': 1141, 'HNRNPA1': 1119, 'CELF2': 1096, 'IGF2BP2': 1011, 'TARDBP': 947, 'DDX3X': 864, 'PTBP1': 726, 'CPSF6': 655, 'ATXN2': 566, 'CSTF2T': 560, 'FIP1L1': 546, 'LIN28B': 529, 'RBP_occupancy': 500, 'FUS': 460, 'CPSF7': 456, 'ACI': 453, 'NUDT21': 347, 'ZC3H7B': 339, 'TIAL1': 335, 'EIF4A3': 252, 'WDR33': 246, 'YTHDC1': 231, 'HNRNPD': 211, 'TDP43': 210, 'EWSR1': 209, 'RBM15': 196, 'TIA1': 190, 'HNRNPU': 189, 'ALYREF': 184, 'SRRM4': 181, 'MOV10': 175, 'TAF15': 167, 'RBM47': 157, 'RBFOX1': 154, 'RBM15B': 147, 'CAPRIN1': 130, 'IGF2BP1': 118, 'QKI': 117, 'SMN': 115, 'KHDRBS1': 113, 'CPSF1': 112, 'YTHDF2': 110, 'IGF2BP3': 104, 'DGCR8': 102, 'FMR1': 98, 'YTHDF1': 90, 'SRSF7': 87, 'SRSF10': 86, 'STAU1': 81, 'SDOS': 76, 'BUD13': 69, 'CPSF3': 66, 'SUB1': 65, 'HNRNPF': 65, 'KHSRP': 64, 'AQR': 64, 'DIS3L2': 63, 'YBX3': 63, 'CPSF4': 60, 'RTCB': 59, 'MBNL2': 58, 'FBL': 57, 'ALKBH1': 56, 'GEMIN5': 55, 'PPIL4': 54, 'PPIG': 51, 'METTL3': 51, 'SAFB2': 51, 'SDAD1': 50, 'ALKBH5': 50, 'NOLC1': 49, 'DDX24': 48, 'HNRNPL': 47, 'LIN28A': 45, 'CHTOP': 45, 'DDX3': 44, 'PRPF8': 42, 'TRA2A': 41, 'EIF3G': 40, 'HNRNPM': 40, 'NONO': 39, 'NCBP2': 39, 'SSB': 38, 'UCHL5': 37, 'AKAP1': 36, 'GRWD1': 35, 'LSM11': 33, 'NCBP3': 33, 'PUM2': 32, 'EIF3D': 32, 'BCLAF1': 32, 'PUM1': 32, 'SLTM': 31, 'PABPC4': 30, 'METAP2': 30, 'SF3B4': 30, 'NOP58': 30, 'EIF3A': 29, 'HNRNPA2B1': 29, 'ZNF800': 28, 'XRN2': 28, 'RBFOX2': 27, 'FTO': 27, 'FXR2': 26, 'SRSF3': 26, 'EFTUD2': 25, 'ZNF622': 24, 'FAM120A': 23, 'DROSHA': 23, 'TARBP2': 23, 'NKRF': 22, 'UTP3': 22, 'YTHDC2': 20, 'CPSF2': 20, 'DDX6': 20, 'AGGF1': 20, 'SF3B1': 20, 'XPO5': 20, 'ZC3H11A': 19, 'EIF3B': 19, 'FASTKD2': 17, 'PABPN1': 17, 'SERBP1': 17, 'MKRN1': 17, 'DHX30': 16, 'LARP4': 16, 'GTF2F1': 15, 'SMNDC1': 14, 'ILF3': 14, 'YTHDF3': 14, 'FUBP3': 14, 'METTL14': 14, 'PCBP2': 14, 'PRPF4': 13, 'YWHAG': 12, 'FXR1': 12, 'ZRANB2': 12, 'TBRG4': 12, 'NOP56': 12, 'AATF': 12, 'GPKOW': 12, 'WRN': 11, 'MTPAP': 11, 'NIPBL': 11, 'LARP7': 11, 'NIP7': 11, 'WTAP': 10, 'PCBP1': 9, 'U2AF1': 9, 'DDX51': 9, 'FKBP4': 9, 'HLTF': 8, 'SUPV3L1': 7, 'EIF4G2': 7, 'HNRNPK': 7, 'SAFB': 7, 'EIF3H': 7, 'SF3A3': 7, 'NOL12': 7, 'AKAP8L': 7, 'RPS3': 6, 'MATR3': 6, 'DDX52': 5, 'NSUN2': 5, 'XRCC6': 5, 'EXOSC5': 5, 'RBM22': 5, 'DKC1': 5, 'DDX42': 4, 'WDR3': 4, 'DDX55': 4, 'TNRC6B': 4, 'TNRC6A': 3, 'RBPMS': 3, 'SFPQ': 3, 'AARS': 3, 'DDX59': 3, 'DDX21': 3, 'SUGP2': 3, 'RBM10': 3, 'TNRC6C': 3, 'UTP18': 2, 'CDC40': 2, 'CPEB4': 2, 'SRSF9': 2, 'GRSF1': 2, 'TROVE2': 2, 'ZC3H8': 2, 'WDR43': 2, 'DICER1': 1, 'UNR': 1, 'ABCF1': 1, 'PHF6': 1, 'BCCIP': 1, 'HNRNPH1': 1, 'RBM5': 1, 'STAU2': 1, 'RPS11': 1}
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
plt.title('Top 20 Protein Counts Triplex Combined')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
