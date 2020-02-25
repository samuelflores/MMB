set(ENV{OPENMM_INCLUDE_PATH} "/home/sam/github/openmm/./include;/home/sam/github/openmm/./include/openmm;/home/sam/github/openmm/./include/openmm/internal;/home/sam/github/openmm/openmmapi/include;/home/sam/github/openmm/openmmapi/include/openmm;/home/sam/github/openmm/openmmapi/include/openmm/internal;/home/sam/github/openmm/olla/include;/home/sam/github/openmm/olla/include/openmm;/home/sam/github/openmm/olla/include/openmm/internal;/home/sam/github/openmm/serialization/include;/home/sam/github/openmm/serialization/include/openmm;/home/sam/github/openmm/serialization/include/openmm/internal;/home/sam/github/openmm/plugins/amoeba/openmmapi/include;/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm;/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal;/home/sam/github/openmm/plugins/rpmd/openmmapi/include;/home/sam/github/openmm/plugins/rpmd/openmmapi/include/openmm;/home/sam/github/openmm/plugins/rpmd/openmmapi/include/openmm/internal;/home/sam/github/openmm/plugins/drude/openmmapi/include;/home/sam/github/openmm/plugins/drude/openmmapi/include/openmm;/home/sam/github/openmm/plugins/drude/openmmapi/include/openmm/internal")
set(ENV{OPENMM_LIB_PATH} "/usr/local/openmm/lib")
message("OPENMM_LIB_PATH = " $ENV{OPENMM_LIB_PATH})
message("OPENMM_INCLUDE_PATH = " $ENV{OPENMM_INCLUDE_PATH})
execute_process(
    COMMAND "/usr/bin/python" setup.py install --root=$ENV{DESTDIR}/
    WORKING_DIRECTORY ""
)
