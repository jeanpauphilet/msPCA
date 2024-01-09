#pragma once
class ConstantArguments
{
public:
    static constexpr int separatorLengthLong = 3;
    static constexpr int wordLengthShort = 10;
    static constexpr int wordLengthMiddle = 20;
    static constexpr int wordLengthLong = 25;
    static constexpr int separatorLengthShort = 1;

    static constexpr int iterationNumberLimitToPrint = 25;
    static constexpr int iterationNumberIntervalToPrint = 10;

    static constexpr int precisionForObjectiveValue = 3;
    static constexpr int precisionForOrthogonalityViolation = 2;
    static constexpr int precisionForTime = 3;

    static constexpr int millisecondsToSeconds = 1000;

    static constexpr int valueInDiagonal = 1;

    static constexpr int defaultMaximumNumberOfIterations = 200;
    static constexpr bool defaultIsVerbose = true;
    static constexpr double defaultViolationTolerance = 1e-4;

    static constexpr double initialOfvBest = -1e10;

    static constexpr double slowPeriodRate = 0.15;
    static constexpr double fastPeriodRate = 0.75;
    static constexpr double changedRateLow = 0.01;
    static constexpr double changedRateHigh = 0.05;

    static constexpr int initialIterationNumber = 1;

    static constexpr int rateForSigmaCurrent = 2;
    static constexpr double differenceForLambda0 = 1e-4;
    static constexpr double lowerLimitForViolation = 1e-7;

    static constexpr int numbersOfSparseVersions = 100;
    static constexpr int appliedTimeLimit = 20;
    static constexpr int defaultTimeLimit = 720;
    static constexpr int defaultSupportValue = 0;
    static constexpr int defaultSupportSize = 1;
    static constexpr int filledSupportValue = -1;
    static constexpr int acceptableSupportValue = 1;
    static constexpr int defaultCountdown = 100;
    static constexpr int differenceFromMimValue = 1;

    static constexpr double ofvPrecision = 1e-8;
};

