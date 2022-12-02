SRC=lb.f90 #${wildcard *.f90}
OBJGPU=${patsubst %.f90,build_gpu/%.o,${SRC}}
OBJCPU=${patsubst %.f90,build_cpu/%.o,${SRC}}
GPU=build_gpu/lbopenacc.exe
CPU=build_cpu/lbopenacc.exe

# nvaccelinfo
GPUFLAGS=-Minfo=accel -ta=tesla:cc86
CPUFLAGS=-Minfo=opt

all: gpu cpu
	@echo "Done."

gpu: ${OBJGPU}
	nvfortran ${OBJGPU} -o ${GPU} ${LDFLAGS} ${GPUFLAGS}

build_gpu/%.o : %.f90
	mkdir -p ${dir $@}
	nvfortran -o $@ $< -c ${GPUFLAGS}

cpu: ${OBJCPU}
	nvfortran ${OBJCPU} -o ${CPU} ${LDFLAGS} ${CPUFLAGS}

build_cpu/%.o : %.f90
	mkdir -p ${dir $@}
	nvfortran -o $@ $< -c ${CPUFLAGS}

clean:
	rm -rf build_gpu/*
	rm -rf build_cpu/*
