import openpyxl as pxl
from win32com.client import Dispatch
import paths, os


#colors
cRed = 		pxl.styles.PatternFill(start_color='FFFF0000', end_color='FFFF0000', fill_type='solid')
cPurple = 	pxl.styles.PatternFill(start_color='FFA313FF', end_color='FFA313FF', fill_type='solid')
cGreen = 	pxl.styles.PatternFill(start_color='FF67FF86', end_color='FF67FF86', fill_type='solid')
cOrange = 	pxl.styles.PatternFill(start_color='FFFFBF65', end_color='FFFFBF65', fill_type='solid')
cGrey = 	pxl.styles.PatternFill(start_color='FFBBBBBB', end_color='FFBBBBBB', fill_type='solid')
status_cols = {'S': cGreen, 'F':cRed, 'R':cPurple, 'W':cOrange, 'C':cGrey}

#border
bSides = pxl.styles.borders.Border(left=pxl.styles.borders.Side(style='thin'), right=pxl.styles.borders.Side(style='thin'))
bHeader = pxl.styles.borders.Border(left=pxl.styles.borders.Side(style='thin'), right=pxl.styles.borders.Side(style='thin'), bottom=pxl.styles.borders.Side(style='double'))

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
		ws.cell(column=col, row=row).font = pxl.styles.Font(bold=True)
		ws.cell(column=col, row=row).border = bHeader
		if n == 'STATUS':
			status_idx = i
		if n == 'DIRECTORY':
			directory_idx = i
		if n == 'OUTXYZ':
			outxyz_idx = i
		if n == 'INXYZ':
			inxyz_idx = i

	directories = [d[directory_idx] for d in data]
	statuses = [d[status_idx] for d in data]
	Nrows = len(data)
	Ncols = len(data[0])
	for i, d in zip(range(Nrows), data):
		row = i + 2
		status_col = status_cols[statuses[i]]
		for j, n, x in zip(range(Ncols), fieldnames, d):
			col = j+1
			cell = ws.cell(column=col, row=row)
			if j == 0:
				cell.font = pxl.styles.Font(bold=True)
			if n == 'DIRECTORY':
				cell.value = f'=HYPERLINK("{os.path.join(paths.master, x)}", "{x}")'
				cell.style = 'Hyperlink'
			elif n == 'OUTXYZ':
				if not x == '':
					cell.value = f'=HYPERLINK("{os.path.join(paths.master,directories[i])}/molview2_start_out.bat", "output.xyz")'
					cell.style = 'Hyperlink'
			elif n == 'INXYZ':
				if not x == '':
					cell.value = f'=HYPERLINK("{os.path.join(paths.master,directories[i])}/molview2_start_in.bat", "input.xyz")'
					cell.style = 'Hyperlink'
			else:
				cell.value = x

			if n in ['ID', 'STATUS']:
				cell.fill = status_col
			cell.border = bSides



last_col = pxl.utils.cell.get_column_letter(col)
ws.auto_filter.ref = f"A1:{last_col}{row}"
ws.freeze_panes = ws['B2']


wb.save(paths.results_table_pretty)

excel = Dispatch('Excel.Application')
wb = excel.Workbooks.Open(paths.results_table_pretty)
excel.Worksheets(1).Activate()
excel.ActiveSheet.Columns.AutoFit()
wb.Save()
wb.Close()