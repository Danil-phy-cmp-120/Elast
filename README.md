# Elast
Code for calculating elastic modulus with VASP

1) Создайте папку "ION". Проведите ионную релаксацию в этой папке используя VASP и сведите внешнее давление к нулю (приемлимое значение порядка 1 кбар). Поместите файлы проекта рядом с папкой "ION". 
2) Запустите программу create_task.py. Данная программа создаст файлы задачи для кода VASP с искаженным базисами в POSCAR.
3) Запстите скрипт vasp_qsub_elast командой "qsub -q node vasp_qsub_elast". Данный скрипт поочередно запустист рассчеты VASP во всех вложенных папках в директории "bulk".
4) После завершения рассчетов запустите программу Read_out.py. Эта программа считает значения полной энергии из файлов OUTCAR и вычислит отдельные коэффициенты тензора упругих постоянных.


