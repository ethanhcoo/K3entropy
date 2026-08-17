
CFLAGS = -g
VERIFY_CFLAGS = -std=c11 -O2 -Wall -Wextra -Wpedantic -Werror
CPPFLAGS = -I/opt/homebrew/include -Icertification -Iillustration_and_exploration
LDFLAGS = -L/opt/homebrew/lib -lmpfr -lgmp
BIN = ../bin
CC = gcc

CERTIFICATION_DIR = certification
ILLUSTRATION_DIR = illustration_and_exploration
CERTIFICATE_COMMON_OBJECTS = \
	${CERTIFICATION_DIR}/certificate_interval.o \
	${CERTIFICATION_DIR}/k3_orbit_data.o
ILLUSTRATION_COMMON_OBJECTS = \
	${ILLUSTRATION_DIR}/k3.o \
	${ILLUSTRATION_DIR}/ps.o

.PHONY: all periodic periodicExact verifyPeriodicOrbit derivatives \
	derivativesExact derivativesInterval arcs mclass verifyArcs \
	testArcWord testCertificateInterval testK3OrbitData \
	checkCertificateCommon checkPathCertification \
	checkDerivativeCertification checkPeriodicOrbitCertification \
	checkAllCertifications clean

all: periodic periodicExact verifyPeriodicOrbit derivatives derivativesExact \
	derivativesInterval arcs mclass verifyArcs

periodic: ${ILLUSTRATION_DIR}/periodic.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} -lm -o ${BIN}/periodic

periodicExact: ${ILLUSTRATION_DIR}/periodicExact.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} -lm -o ${BIN}/periodicExact

# Outward-rounded replacement for periodicExact's pseudo-orbit check.
verifyPeriodicOrbit: ${CERTIFICATION_DIR}/verifyPeriodicOrbit.o ${CERTIFICATE_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} -lm -o ${BIN}/verifyPeriodicOrbit

# Recompute all domain, branch, residual, and Appendix C.3 table checks.
checkPeriodicOrbitCertification: verifyPeriodicOrbit
	${BIN}/verifyPeriodicOrbit

derivatives: ${ILLUSTRATION_DIR}/derivatives.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} -lm -o ${BIN}/derivatives

derivativesExact: ${ILLUSTRATION_DIR}/derivativesExact.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} -lm -o ${BIN}/derivativesExact

derivativesInterval: ${CERTIFICATION_DIR}/derivativesInterval.o ${CERTIFICATE_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} -lm -o ${BIN}/derivativesInterval

# Run every certificate described in Appendix C.2; any inconclusive check
# makes derivativesInterval exit unsuccessfully.
checkDerivativeCertification: derivativesInterval
	${BIN}/derivativesInterval

arcs: ${ILLUSTRATION_DIR}/arcs.o ${ILLUSTRATION_COMMON_OBJECTS}
	${CC} $^ -lm -o ${BIN}/arcs

mclass: ${ILLUSTRATION_DIR}/mclass.o ${ILLUSTRATION_COMMON_OBJECTS} \
	${CERTIFICATION_DIR}/arc_certificate_data.o
	${CC} $^ -lm -o ${BIN}/mclass

verifyArcs: ${CERTIFICATION_DIR}/verifyArcs.o \
	${CERTIFICATION_DIR}/arc_word.o \
	${CERTIFICATION_DIR}/arc_certificate_data.o \
	${CERTIFICATE_COMMON_OBJECTS}
	${CC} $^ ${LDFLAGS} -lm -o ${BIN}/verifyArcs

testArcWord: ${CERTIFICATION_DIR}/testArcWord.o \
	${CERTIFICATION_DIR}/arc_word.o \
	${CERTIFICATION_DIR}/arc_certificate_data.o
	${CC} $^ -o ${BIN}/testArcWord

testCertificateInterval: ${CERTIFICATION_DIR}/testCertificateInterval.o \
	${CERTIFICATION_DIR}/certificate_interval.o
	${CC} $^ ${LDFLAGS} -o ${BIN}/testCertificateInterval

testK3OrbitData: ${CERTIFICATION_DIR}/testK3OrbitData.o \
	${CERTIFICATION_DIR}/k3_orbit_data.o
	${CC} $^ -o ${BIN}/testK3OrbitData

checkCertificateCommon: testCertificateInterval testK3OrbitData
	${BIN}/testCertificateInterval
	${BIN}/testK3OrbitData

checkPathCertification: verifyArcs testArcWord
	${BIN}/testArcWord
	${BIN}/verifyArcs

checkAllCertifications: checkCertificateCommon \
	checkPeriodicOrbitCertification checkDerivativeCertification \
	checkPathCertification

# Certification sources always compile with the strict certificate flags.
${CERTIFICATION_DIR}/%.o: ${CERTIFICATION_DIR}/%.c
	${CC} ${CPPFLAGS} ${VERIFY_CFLAGS} -c $< -o $@

# Illustration and exploratory programs retain their original debug build.
${ILLUSTRATION_DIR}/%.o: ${ILLUSTRATION_DIR}/%.c
	${CC} ${CPPFLAGS} ${CFLAGS} -c $< -o $@

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
