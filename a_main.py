# -*- coding: utf-8 -*-
import sys
import os

# Pandas
import pandas as pd

# Partial functool
from functools import partial

# Pyqt5
from PyQt5 import QtCore, QtGui, QtWidgets

#Lug Class
import b_lug_class as lc

# Gui Script
from c_gui import Ui_MainWindow


class main_ui(Ui_MainWindow): 


	def __init__(self,MainWindow): 

		Ui_MainWindow.setupUi(self,MainWindow)

		
		# Bearing tab
		# ---------------------------------------------------------------------------------------
		# read  excel as df
		self.bearing_toolButton_excel_input.clicked.connect(partial(self.set_input_excel,"bearing"))
		
		# Calulate 
		self.bearing_pushButton.clicked.connect(partial(self.calculate,"bearing"))


		# Plain Buching
		# ---------------------------------------------------------------------------------------
		# read excel as df
		self.plain_bushing_toolButton_excel_input.clicked.connect(partial(self.set_input_excel,"plain_bushing"))
		
		# Calulate 
		self.plain_bushing_pushButton.clicked.connect(partial(self.calculate,"plain_bushing"))


		# Flange Buching
		# ---------------------------------------------------------------------------------------
		# read excel as df
		self.flange_bushing_toolButton_excel_input.clicked.connect(partial(self.set_input_excel,"flange_bushing"))
		
		# Calulate 
		self.flange_bushing_pushButton.clicked.connect(partial(self.calculate,"flange_bushing"))


		# Pin Bolt
		# ---------------------------------------------------------------------------------------
		# read lug excel as df
		self.pin_bolt_toolButton_excel_input.clicked.connect(partial(self.set_input_excel,"pin_bolt"))
		
		# Calulate 
		self.pin_bolt_pushButton.clicked.connect(partial(self.calculate,"pin_bolt"))
		

		# Lug tab
		# ---------------------------------------------------------------------------------------

		# read lug excel as df
		self.lug_toolButton_excel_input.clicked.connect(partial(self.set_input_excel,"lug"))

		# Calulate 
		self.lug_pushButton.clicked.connect(partial(self.calculate,"lug"))



	## Calculate function
	def calculate(self,element_type):

		# Converting to dict transposing it
		self.elements_dict = self.df.T.to_dict()

		# Poulating lug field with the last item on the excel
		for element in self.elements_dict.keys():

			for p in self.elements_dict[element].keys():
				# Writing input data into cells
				getattr(self, element_type + "_lineEdit_" + p).setText(str(self.elements_dict[element][p]))


		# iterating over the dictionary
		for element in self.elements_dict.keys():

			print(element)
			print("=======================================")
			print("")

			# Initialize Lug class to allocate attributes
			l = getattr(lc, element_type)()

			for p in self.elements_dict[element].keys():

				# Writing input data into cells
				getattr(self, element_type + "_lineEdit_" + p).setText(str(self.elements_dict[element][p]))

				# Reading modified lugs dict
				self.elements_dict[element][p] = float(str(getattr(self, element_type + "_lineEdit_" + p).text()))

				# Modifiying attributes of the class
				setattr(l, p, self.elements_dict[element][p])
				
			# Calculate Lug
			l.calculate()

			print("")
			print("")

			for p in l.__dict__.keys():
				self.elements_dict[element][p] = round(l.__dict__[p],2)


		self.set_results_table(self.elements_dict, table=getattr(self, element_type + "_tableWidget_results")) 
		self.write_to_excel(self.elements_dict,element_type)


	# GUI Handle Functions 
	# ============================================================================================

	def set_input_excel(self,element):


		# Excel path
		self.dlg = QtWidgets.QFileDialog()
		getattr(self, element + "_lineEdit_excel_input").setText(str(self.dlg.getOpenFileName()[0]))
		self.path_excel = str(getattr(self, element + "_lineEdit_excel_input").text())


		# Read Excel File
		self.df = pd.read_excel(self.path_excel, sheet_name=element, index_col = 0)

		print(self.df)



	def set_results_table(self,sol_dict,table):

		# Allocating the number of results in the column
		table.setRowCount(len(sol_dict))
		
		# Renaming columns of the gui table
		table.setColumnCount(len(sol_dict[list(sol_dict.keys())[0]]))
		table.setHorizontalHeaderLabels(sol_dict[list(sol_dict.keys())[0]].keys())


		i = 0
		for m in sol_dict.keys():
			j = 0
			for p in sol_dict[m]: 
				table.setItem(i, j, QtWidgets.QTableWidgetItem(str(sol_dict[m][p])))
				j += 1
			i += 1

	def write_to_excel(self,sol_dict,element):

		# Load package to load workbooks
		from openpyxl import load_workbook

		# Coverting dict to Dataframe
		df = pd.DataFrame.from_dict(sol_dict)

		book = load_workbook("Output"+".xlsx")
		writer = pd.ExcelWriter("Output" +".xlsx", engine = "openpyxl")
		writer.book = book 
		writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
		df.T.to_excel(writer,element)
		writer.save()



if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = main_ui(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
