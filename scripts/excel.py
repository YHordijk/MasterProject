import openpyxl as pxl
from win32com.client import Dispatch



wb = pxl.Workbook()
wb.remove(wb.active) #remove default sheet
wb.create_sheet(title='Overview')
ws = wb.get_sheet_by_name('Overview')
with open(r"D:\Users\Yuman\Desktop\MasterProject\results\master_results_table.csv") as csv:
	lines = csv.readlines()
	fieldnames = [n.strip() for n in lines[1].split(',')]
	data = [[x.strip() for x in d.split(',')] for d in lines[2:]]

	for i, n in enumerate(fieldnames):
		row = 1
		col = i+1
		ws.cell(column=col, row=row, value=n)

	for i, d in enumerate(data):
		row = i + 2
		for j, x in enumerate(d):
			col = j+1
			ws.cell(column=col, row=row, value=x)



wb.save(r"D:\Users\Yuman\Desktop\MasterProject\results\results_formatted.xlsx")

excel = Dispatch('Excel.Application')
wb = excel.Workbooks.Open(r"D:\Users\Yuman\Desktop\MasterProject\results\results_formatted.xlsx")
excel.Worksheets(1).Activate()
excel.ActiveSheet.Columns.AutoFit()
wb.Save()
wb.Close()