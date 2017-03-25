function lab1Main(dataFile)


arrayOfLines = fileImport(dataFile);


[numberOfForces, forcesCoordinates, externalForces, externalMoments, locationSupports, uVecReactionDir, momentOrForce] = dataManipulation(arrayOfLines);


[A,b] = aBMatrixManipulation(numberOfForces, forcesCoordinates, externalForces, externalMoments, locationSupports, uVecReactionDir, momentOrForce);


x = A\b




















