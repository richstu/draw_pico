#!/bin/env python3
import os

Import("exportEnv")
envClone = exportEnv.Clone()

# source_directories = ['core', 'higgsino', ...]
source_directories = next(os.walk('src'))[1]
if len(source_directories) == 0:
  source_directories = ['./']

envClone.Append(CPPPATH = ['#/inc'])

# Determine names for trees
tree_variable_filenames = [os.path.basename(str(x)) for x in Glob('#/txt/variables/*')]
tree_variable_files = []
tree_generated_files = ['src/core/baby.cpp','inc/core/baby.hpp']
tree_generated_files_mark = []
for tree_name in tree_variable_filenames: # tree_name = 'pico'
  tree_variable_files.append('txt/variables/'+tree_name)
  tree_generated_files.extend([
                          'src/core/baby_'+tree_name+'.cpp',
                          'inc/core/baby_'+tree_name+'.hpp',
                         ])
tree_generated_files_mark = ['#/'+term for term in tree_generated_files]

# Make tree generator
tree_generator_file = 'src/core/generate_baby.cxx'
tree_generator_exe = 'kernel/'+envClone['kernel']+'/run/core/generate_baby.exe'
tree_generator = envClone.Program('#/'+tree_generator_exe, Glob(tree_generator_file))

# Run tree generator
run_tree_generator = envClone.Command(tree_generated_files_mark, [], './'+tree_generator_exe+' '+' '.join(tree_variable_filenames))
envClone.Depends(run_tree_generator, tree_generator_file)
envClone.Depends(run_tree_generator, tree_variable_files)
envClone.Depends(run_tree_generator, tree_generator)

# Run genreflex
reflex_files = ['RooMultiPdf']
reflex_generated_files = []
run_genreflex = []
for reflex_file in reflex_files:
  in_filename = 'core/'+reflex_file+'.hpp'
  out_filenames = ['src/core/'+reflex_file+'_dict.cpp',
                   'src/core/'+reflex_file+'_dict_rdict.pcm']
  reflex_generated_files += out_filenames
  run_genreflex.append(envClone.Command(out_filenames, [], 'cd inc/ && genreflex '+in_filename+' -o ../$TARGET && cd -'))
  envClone.Depends(run_genreflex[-1], 'inc/'+in_filename)
  #	genreflex $(SRC_DIR)/classes.h -s $(SRC_DIR)/classes_def.xml -o $(OBJ_DIR)/a/$(DICTNAME).cc --deep --fail_on_warnings --rootmap=$(OBJ_DIR)/a/$(DICTNAME).rootmap --rootmap-lib=$(SONAME) -I$(PARENT_DIR) 

# Make libraries
libraries = {}
envClone.Append(LINKFLAGS='-L'+'lib/'+envClone['kernel'])

# Make core library
core_objects = []
# Build baby objects
for source_file in tree_generated_files:
  if 'hpp' in source_file: continue
  source_object = envClone.SharedObject(source_file.replace('cpp','os'), '#/'+source_file)
  core_objects.append(source_object)
# Build reflex objects
for source_file in reflex_generated_files:
  if 'pcm' in source_file: continue
  reflex_flags = '-O2 -isystem `root-config --incdir` `root-config --cflags`'
  source_object = envClone.SharedObject(source_file.replace('cpp','os'), source_file, CCFLAGS=reflex_flags)
  core_objects.append(source_object)
# Build non-baby non-reflex core objects
for lib_file in Glob("src/core/*.cpp", exclude=tree_generated_files+reflex_generated_files): 
  source_object = envClone.SharedObject(lib_file)
  core_objects.append(source_object)
# Make core library
libraryName = "DrawPicoCore"
library = envClone.SharedLibrary(target='#/lib/'+envClone['kernel']+'/'+libraryName, source=core_objects)
libraries[libraryName]= library

# Make other libraries
for source_directory in source_directories:
  if source_directory == "core": continue # already made a library for core
  if source_directory == "pythonbindings": continue # do this separately
  libraryName = "DrawPico"+source_directory.capitalize()
  source_files = set()
  for lib_file in Glob("src/"+source_directory+"/*.cpp"): source_files.add(lib_file)
  library = envClone.SharedLibrary(target='#/lib/'+envClone['kernel']+'/'+libraryName, source=sorted(source_files), LIBS=["DrawPicoCore"])
  envClone.Depends(library, libraries["DrawPicoCore"])
  libraries[libraryName]= library

#Make python bindings library 
#this depends on all the others so must be done last
libraryName = "DrawPicoPythonbindings"
source_files = set()
for lib_file in Glob("src/pythonbindings/*.cpp"): source_files.add(lib_file)
allLibraryNames = []
for source_directory in source_directories:
  if source_directory == "pythonbindings": continue
  allLibraryNames.append("DrawPico"+source_directory.capitalize())
library = envClone.SharedLibrary(target='#/lib/'+envClone['kernel']+'/'+libraryName, source=sorted(source_files), LIBS=allLibraryNames)
for libraryName in allLibraryNames:
  envClone.Depends(library, libraries[libraryName])
libraries[libraryName]= library

# Make binaries for every directory
exclude_files = [tree_generator_file] + tree_generated_files
for source_directory in source_directories: # source_directories = 'core'
  for source_file in Glob("src/"+source_directory+"/*.cxx", exclude=exclude_files):
    # Build object
    source_object = envClone.Object(source_file)
    if source_directory!='core': envClone.Depends(source_object, libraries['DrawPicoCore'])
    libraryName = "DrawPico"+source_directory.capitalize()
    # Make program
    program = envClone.Program('#/kernel/'+envClone['kernel']+'/run/'+source_directory+'/${SOURCE.filebase}.exe', source_object, LIBS=["DrawPicoCore", libraryName])
    envClone.Depends(program, libraries["DrawPicoCore"])
    envClone.Depends(program, libraries[libraryName])
    # Make script that links to binary
    source_basename = os.path.splitext(os.path.basename(str(source_file)))[0]
    source_script_path = 'run/'+source_directory+'/'+source_basename+'.exe'
    source_script_path_mark = "#/"+source_script_path
    envClone.Command(source_script_path_mark, program, './scripts/make_run_scripts.py '+source_script_path)
