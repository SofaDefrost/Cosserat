#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Younes'
__copyright__ = '(c) {year}, Inria, {project_name}'
__credits__ = ['{credit_list}']
__license__ = '{license}'
__version__ = '{mayor}.{minor}.{rel}'
__maintainer__ = '{maintainer}'
__email__ = 'yinoussa.adagolodjo-AT-inria.fr'
__status__ = '{dev_status}'
__date__ = "{Date}"

"""
{Description}
{License_info}
"""
import PySimpleGUI as sg

sg.theme('DarkAmber')  # Add a touch of color
# All the stuff inside your window.
layout = [[sg.Text('Some text on Row 1')],
          [sg.Text('Enter something on Row 2'), sg.InputText()],
          [sg.Button('Ok'), sg.Button('Cancel')]]

# Create the Window
window = sg.Window('Window Title', layout)
# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Cancel':  # if user closes window or clicks cancel
        break
    print('You entered ', values[0])

window.close()
