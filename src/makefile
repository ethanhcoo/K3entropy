
CC ?= cc
CFLAGS ?= -g
VERIFY_CFLAGS = -std=c11 -O2 -Wall -Wextra -Wpedantic -Werror
PROJECT_CPPFLAGS = -Icertification -Iillustration_and_exploration
LDLIBS ?= -lmpfr -lgmp
BIN ?= ../bin
PYTHON ?= python3
MAPPING_CLASS_WORD ?= ../mapping_class_word.txt
MANUSCRIPT ?= ../../main_updated.tex

CERTIFICATION_DIR = certification
ILLUSTRATION_DIR = illustration_and_exploration
CERTIFICATE_COMMON_OBJECTS = \
	${CERTIFICATION_DIR}/certificate_interval.o \
	${CERTIFICATION_DIR}/k3_orbit_data.o
ILLUSTRATION_COMMON_OBJECTS = \
	${ILLUSTRATION_DIR}/k3.o \
	${ILLUSTRATION_DIR}/ps.o
PROGRAMS = periodic periodicExact verifyPeriodicOrbit derivatives \
	derivativesExact derivativesInterval arcs mclass verifyArcs \
	testArcWord testCertificateInterval testK3OrbitData

.PHONY: all ${PROGRAMS} \
	checkCertificateCommon checkPathCertification \
	checkDerivativeCertification checkPeriodicOrbitCertification \
	checkAllCertifications checkMappingClassWord checkManuscriptWord clean

all: periodic periodicExact verifyPeriodicOrbit derivatives derivativesExact \
	derivativesInterval arcs mclass verifyArcs

${PROGRAMS}: | ${BIN}

${BIN}:
	mkdir -p $@

periodic: ${ILLUSTRATION_DIR}/periodic.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} ${LDLIBS} -lm -o ${BIN}/periodic

periodicExact: ${ILLUSTRATION_DIR}/periodicExact.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} ${LDLIBS} -lm -o ${BIN}/periodicExact

# Outward-rounded replacement for periodicExact's pseudo-orbit check.
verifyPeriodicOrbit: ${CERTIFICATION_DIR}/verifyPeriodicOrbit.o ${CERTIFICATE_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} ${LDLIBS} -lm -o ${BIN}/verifyPeriodicOrbit

# Recompute all domain, branch, residual, and Appendix C.3 table checks.
checkPeriodicOrbitCertification: verifyPeriodicOrbit
	${BIN}/verifyPeriodicOrbit

derivatives: ${ILLUSTRATION_DIR}/derivatives.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} ${LDLIBS} -lm -o ${BIN}/derivatives

derivativesExact: ${ILLUSTRATION_DIR}/derivativesExact.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} ${LDLIBS} -lm -o ${BIN}/derivativesExact

derivativesInterval: ${CERTIFICATION_DIR}/derivativesInterval.o ${CERTIFICATE_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} ${LDLIBS} -lm -o ${BIN}/derivativesInterval

# Run every certificate described in Appendix C.2; any inconclusive check
# makes derivativesInterval exit unsuccessfully.
checkDerivativeCertification: derivativesInterval
	${BIN}/derivativesInterval

arcs: ${ILLUSTRATION_DIR}/arcs.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} -lm -o ${BIN}/arcs

mclass: ${ILLUSTRATION_DIR}/mclass.o ${ILLUSTRATION_COMMON_OBJECTS} \
	${CERTIFICATION_DIR}/arc_certificate_data.o
	${CC} $^ ${LDFLAGS} -lm -o ${BIN}/mclass

verifyArcs: ${CERTIFICATION_DIR}/verifyArcs.o \
	${CERTIFICATION_DIR}/arc_word.o \
	${CERTIFICATION_DIR}/arc_certificate_data.o \
	${CERTIFICATE_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} ${LDLIBS} -lm -o ${BIN}/verifyArcs

testArcWord: ${CERTIFICATION_DIR}/testArcWord.o \
	${CERTIFICATION_DIR}/arc_word.o \
	${CERTIFICATION_DIR}/arc_certificate_data.o
	${CC} $^ ${LDFLAGS} -o ${BIN}/testArcWord

testCertificateInterval: ${CERTIFICATION_DIR}/testCertificateInterval.o \
	${CERTIFICATION_DIR}/certificate_interval.o
	${CC} $^ ${LDFLAGS} ${LDLIBS} -o ${BIN}/testCertificateInterval

testK3OrbitData: ${CERTIFICATION_DIR}/testK3OrbitData.o \
	${CERTIFICATION_DIR}/k3_orbit_data.o
	${CC} $^ ${LDFLAGS} -o ${BIN}/testK3OrbitData

checkCertificateCommon: testCertificateInterval testK3OrbitData
	${BIN}/testCertificateInterval
	${BIN}/testK3OrbitData

checkPathCertification: verifyArcs testArcWord
	${BIN}/testArcWord
	${BIN}/verifyArcs

checkAllCertifications: checkCertificateCommon \
	checkPeriodicOrbitCertification checkDerivativeCertification \
	checkPathCertification

# This downstream reproduction check is intentionally separate from the
# rigorous C certification aggregate above.
checkMappingClassWord: mclass ${MAPPING_CLASS_WORD}
	${PYTHON} ../check_mapping_class_word.py \
		--word "${MAPPING_CLASS_WORD}" \
		--mclass "${BIN}/mclass"

# Optional cross-check against the manuscript equation carrying the stable
# eq:f2-flipper-word label. Override MANUSCRIPT when it lives elsewhere.
checkManuscriptWord: checkMappingClassWord
	${PYTHON} ../check_mapping_class_word.py \
		--word "${MAPPING_CLASS_WORD}" \
		--manuscript "${MANUSCRIPT}"

# Certification sources always compile with the strict certificate flags.
${CERTIFICATION_DIR}/%.o: ${CERTIFICATION_DIR}/%.c
	${CC} ${PROJECT_CPPFLAGS} ${CPPFLAGS} ${VERIFY_CFLAGS} -c $< -o $@

# Illustration and exploratory programs retain their original debug build.
${ILLUSTRATION_DIR}/%.o: ${ILLUSTRATION_DIR}/%.c
	${CC} ${PROJECT_CPPFLAGS} ${CPPFLAGS} ${CFLAGS} -c $< -o $@

${CERTIFICATION_DIR}/verifyPeriodicOrbit.o: \
	${CERTIFICATION_DIR}/certificate_interval.h \
	${CERTIFICATION_DIR}/k3_orbit_data.h
${CERTIFICATION_DIR}/derivativesInterval.o: \
	${CERTIFICATION_DIR}/certificate_interval.h \
	${CERTIFICATION_DIR}/k3_orbit_data.h
${CERTIFICATION_DIR}/verifyArcs.o: \
	${CERTIFICATION_DIR}/arc_word.h \
	${CERTIFICATION_DIR}/arc_certificate_data.h \
	${CERTIFICATION_DIR}/certificate_interval.h \
	${CERTIFICATION_DIR}/k3_orbit_data.h
${CERTIFICATION_DIR}/certificate_interval.o: ${CERTIFICATION_DIR}/certificate_interval.h
${CERTIFICATION_DIR}/k3_orbit_data.o: ${CERTIFICATION_DIR}/k3_orbit_data.h
${CERTIFICATION_DIR}/arc_word.o: ${CERTIFICATION_DIR}/arc_word.h
${CERTIFICATION_DIR}/arc_certificate_data.o: ${CERTIFICATION_DIR}/arc_certificate_data.h
${CERTIFICATION_DIR}/testArcWord.o: \
	${CERTIFICATION_DIR}/arc_word.h \
	${CERTIFICATION_DIR}/arc_certificate_data.h
${CERTIFICATION_DIR}/testCertificateInterval.o: \
	${CERTIFICATION_DIR}/certificate_interval.h
${CERTIFICATION_DIR}/testK3OrbitData.o: ${CERTIFICATION_DIR}/k3_orbit_data.h

${ILLUSTRATION_DIR}/k3.o: ${ILLUSTRATION_DIR}/k3.h
${ILLUSTRATION_DIR}/ps.o: ${ILLUSTRATION_DIR}/k3.h
${ILLUSTRATION_DIR}/periodic.o: ${ILLUSTRATION_DIR}/k3.h
${ILLUSTRATION_DIR}/periodicExact.o: ${ILLUSTRATION_DIR}/k3.h
${ILLUSTRATION_DIR}/derivatives.o: ${ILLUSTRATION_DIR}/k3.h
${ILLUSTRATION_DIR}/derivativesExact.o: ${ILLUSTRATION_DIR}/k3.h
${ILLUSTRATION_DIR}/arcs.o: ${ILLUSTRATION_DIR}/k3.h
${ILLUSTRATION_DIR}/mclass.o: \
	${ILLUSTRATION_DIR}/k3.h \
	${CERTIFICATION_DIR}/arc_certificate_data.h

clean:
	rm -f ${CERTIFICATION_DIR}/*.o ${ILLUSTRATION_DIR}/*.o *.o ${BIN}/*
