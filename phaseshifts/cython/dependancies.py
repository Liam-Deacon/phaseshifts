from modulefinder import ModuleFinder

finder = ModuleFinder()
finder.run_script('phsh.pyx')
print finder.report()
