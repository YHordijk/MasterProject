import openpyxl as pxl
from win32com.client import Dispatch
import paths


#colors
cRed = 		pxl.styles.PatternFill(start_color='FFFF0000', end_color='FFFF0000', fill_type='solid')
cPurple = 	pxl.styles.PatternFill(start_color='FFA313FF', end_color='FFA313FF', fill_type='solid')
cGreen = 	pxl.styles.PatternFill(start_color='FFE8FFDF', end_color='FFE8FFDF', fill_type='solid')
cOrange = 	pxl.styles.PatternFill(start_color='FFFFAD00', end_color='FFFFAD00', fill_type='solid')
cGrey = 	pxl.styles.PatternFill(start_color='FF898989', end_color='FF898989', fill_type='solid')

status_cols = {'S': cGreen, 'F':cRed, 'R':cPurple, 'W':cOrange, 'C':cGrey}

wb = pxl.Workbook()
wb.remove(wb.active) #remove default sheet
wb.create_sheet(title='Overview')
ws = wb['Overview']
with open(paths.results_table) as csv:
	lines = csv.readlines()
	fieldnames = [n.strip() for n in lines[1].split(',')]
	data = [[x.strip() for x in d.split(',')] for d in lines[2:]]

	for i, n in enumerate(fieldnames):
		row = 1
		col = i+1
		ws.cell(column=col, row=row, value=n)
		if n == 'STATUS':
			status_idx = i

	statuses = [d[status_idx] for d in data]
	for i, d in enumerate(data):
		row = i + 2
		status_col = status_cols[statuses[i]]
		for j, x in enumerate(d):
			col = j+1
			ws.cell(column=col, row=row, value=x)
			ws.cell(column=col, row=row).fill = status_col



wb.save(paths.results_table_pretty)

excel = Dispatch('Excel.Application')
wb = excel.Workbooks.Open(paths.results_table_pretty)
excel.Worksheets(1).Activate()
excel.ActiveSheet.Columns.AutoFit()
wb.Save()
wb.Close()