import xlsxwriter as xw

def writeXLSXProtein(workbook, protein_name, bench_data, atom_data, centroid_data):
    worksheet = workbook.add_worksheet(protein_name)
    bold = workbook.add_format({"bold": 1})
    headings = ["Clique Sizes", "Bench", "New"]
    worksheet.write_row("A1", headings, bold)
    worksheet.write_column("A2", [1,2,3,4,5,6,7])
    worksheet.write_column("B2", bench_data)
    worksheet.write_column("C2", atom_data)
    worksheet.write_column("D2", centroid_data)
    chart = workbook.add_chart({"type": "column"})
    chart.add_series(
    {
        "name": "="+protein_name+"!$B$1",
        "categories": "="+protein_name+"!$A$2:$A$8",
        "values": "="+protein_name+"!$B$2:$B$8"
    })
    chart.add_series(
    {
        "name": "="+protein_name+"!$C$1",
        "categories": "="+protein_name+"!$A$2:$A$8",
        "values": "="+protein_name+"!$C$2:$C$8"
    })
    chart.add_series(
    {
        "name": "="+protein_name+"!$D$1",
        "categories": "="+protein_name+"!$A$2:$A$8",
        "values": "="+protein_name+"!$D$2:$D$8"
    })
    chart.set_style(21)
    worksheet.insert_chart("F2", chart, {"x_offset": 50, "y_offset": 35, "x_scale": 2.2, "y_scale": 1.8})
    worksheet.set_first_sheet()

def writeXLSXCode(workbook, protein_name, codes_bench, codes_atom):
    worksheet = workbook.add_worksheet(protein_name)
    bold = workbook.add_format({"bold": 1})
    headings = ["Bench", "Atom"]
    worksheet.write_row("A1", headings, bold)
    worksheet.write_column("A2", codes_bench)
    worksheet.write_column("B2", codes_atom)
    chart = workbook.add_chart({"type": "column"})
    chart.add_series(
    {
        "name": "="+protein_name+"!$A$1",
        "categories": "="+protein_name+"!$A$1",
        "values": "="+protein_name+"!$A$2"
    })
    chart.add_series(
    {
        "name": "="+protein_name+"!$B$1",
        "categories": "="+protein_name+"!$B$1",
        "values": "="+protein_name+"!$B$2"
    })
    chart.set_style(21)
    worksheet.insert_chart("F2", chart, {"x_offset": 0, "y_offset": 0, "x_scale": 1.5, "y_scale": 1})
    worksheet.set_first_sheet()


    
