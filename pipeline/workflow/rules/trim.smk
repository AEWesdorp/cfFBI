# unit name = unique name?
units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "library_prep": str})
    .set_index(["sample_name"], drop=False)
    .sort_index()
)

SRSLY_index = (
    pd.read_csv("resources/adapter_indexes/SRSLY_index.txt", sep="\t", dtype={"name": str})
    .set_index(["name"], drop=False)
    .sort_index()
)

IDT384UMI_index = (
    pd.read_csv("resources/adapter_indexes/IDT384UMI_index.txt", sep="\t", dtype={"name": str})
    .set_index(["name"], drop=False)
    .sort_index()
)

KAPA_index = (
    pd.read_csv("resources/adapter_indexes/KAPA_index.txt", sep="\t", dtype={"name": str})
    .set_index(["name"], drop=False)
    .sort_index()
)

def get_adapter_R1(sample_name):
    adapter_type = units.loc[sample_name]['adapter_type'].values[0]
    if adapter_type == 'SRSLY_dual_index':
        library_prep = units.loc[sample_name]['library_prep'].values[0]
        UDI = units.loc[sample_name]['UDI']
        if library_prep == "SRSLY":
            adapter_index_P5_R2 = SRSLY_index.loc[UDI]['srsly_i5_rev'].values[0]
            adapter_R2 = config['SRSLY_dual_index']['adapterP5_part1'] + adapter_index_P5_R2 + config['SRSLY_dual_index']['adapterP5_part2']

            adapter_index_P7_R1 = SRSLY_index.loc[UDI]['srsly_i7_rev'].values[0]
            adapter_R1 = config['SRSLY_dual_index']['adapterP7_part1'] + adapter_index_P7_R1 + config['SRSLY_dual_index']['adapterP7_part2']
        elif library_prep != "SRSLY":
            print("ERROR: ibrary_prep is not in line with adapter_type")        

    elif adapter_type == 'IDT384UMI_dual':
        library_prep = units.loc[sample_name]['library_prep'].values[0]
        UDI = units.loc[sample_name]['UDI']
        if library_prep == "KAPA":
            adapter_index_P5_R2 = IDT384UMI_index.loc[UDI]['IDT384UMI_i5_rev'].values[0]
            adapter_R2 = config['IDT384UMI_dual_index']['adapterP5_part1'] + adapter_index_P5_R2 + config['IDT384UMI_dual_index']['adapterP5_part2']

            adapter_index_P7_R1 = IDT384UMI_index.loc[UDI]['IDT384UMI_i7_rev'].values[0]
            adapter_R1 = config['IDT384UMI_dual_index']['adapterP7_part1'] + adapter_index_P7_R1 + config['IDT384UMI_dual_index']['adapterP7_part2']
        elif library_prep != "KAPA":
            print("ERROR: ibrary_prep is not in line with adapter_type") 

    elif adapter_type == 'KAPA_single_index':
        library_prep = units.loc[sample_name]['library_prep'].values[0]
        UDI = units.loc[sample_name]['UDI']
        if library_prep == "KAPA":
            adapter_R2 = config['KAPA_single_index']['adapterP5']

            adapter_index_P7_R1 = KAPA_index.loc[UDI]['kapa_i7'].values[0]
            adapter_R1 = config['KAPA_single_index']['adapterP7_part1'] + adapter_index_P7_R1 + config['KAPA_single_index']['adapterP7_part2']
        elif library_prep != "KAPA":
            print("ERROR: library_prep is not in line with adapter_type")
    else:
        print(f"ERROR: Sample: {sample_name} : adapter_type is not recognized. Adapter type: {adapter_type} \nCan be either 'SRSLY_dual_index' or 'KAPA_single_index'")

    return adapter_R1


def get_adapter_R2(sample_name):
    adapter_type = units.loc[sample_name]['adapter_type'].values[0]
    if adapter_type == 'SRSLY_dual_index':
        library_prep = units.loc[sample_name]['library_prep'].values[0]
        UDI = units.loc[sample_name]['UDI']
        if library_prep == "SRSLY":
            adapter_index_P5_R2 = SRSLY_index.loc[UDI]['srsly_i5_rev'].values[0]
            adapter_R2 = config['SRSLY_dual_index']['adapterP5_part1'] + adapter_index_P5_R2 + config['SRSLY_dual_index']['adapterP5_part2']

            adapter_index_P7_R1 = SRSLY_index.loc[UDI]['srsly_i7_rev'].values[0]
            adapter_R1 = config['SRSLY_dual_index']['adapterP7_part1'] + adapter_index_P7_R1 + config['SRSLY_dual_index']['adapterP7_part2']
        elif library_prep != "SRSLY":
            print("ERROR: library_prep is not in line with adapter_type")    


    elif adapter_type == 'IDT384UMI_dual':
        library_prep = units.loc[sample_name]['library_prep'].values[0]
        UDI = units.loc[sample_name]['UDI']
        if library_prep == "KAPA":
            adapter_index_P5_R2 = IDT384UMI_index.loc[UDI]['IDT384UMI_i5_rev'].values[0]
            adapter_R2 = config['IDT384UMI_dual_index']['adapterP5_part1'] + adapter_index_P5_R2 + config['IDT384UMI_dual_index']['adapterP5_part2']

            adapter_index_P7_R1 = IDT384UMI_index.loc[UDI]['IDT384UMI_i7_rev'].values[0]
            adapter_R1 = config['IDT384UMI_dual_index']['adapterP7_part1'] + adapter_index_P7_R1 + config['IDT384UMI_dual_index']['adapterP7_part2']
        elif library_prep != "KAPA":
            print("ERROR: ibrary_prep is not in line with adapter_type")    

    elif adapter_type == 'KAPA_single_index':
        library_prep = units.loc[sample_name]['library_prep'].values[0]
        UDI = units.loc[sample_name]['UDI']
        if library_prep == "KAPA":
            adapter_R2 = config['KAPA_single_index']['adapterP5']

            adapter_index_P7_R1 = KAPA_index.loc[UDI]['kapa_i7'].values[0]
            adapter_R1 = config['KAPA_single_index']['adapterP7_part1'] + adapter_index_P7_R1 + config['KAPA_single_index']['adapterP7_part2']
        elif library_prep != "KAPA":
            print("ERROR: library_prep is not in line with adapter_type")
    else:
        print(f"ERROR: Sample: {sample_name} : adapter_type is not recognized. Adapter type: {adapter_type} \nCan be either 'SRSLY_dual_index' or 'KAPA_single_index'")
    return adapter_R2
