// Utility functions
function reverseStr(str) {
    return str.split('').reverse().join('');
}

function countBase(seq, base) {
    return (seq.match(new RegExp(base, 'g')) || []).length;
}

function loopDetective(duplexStr) {
    let [seq1, seq2] = duplexStr.split("\n");
    let matchAscLi = [149, 138];
    let gapAscLi = [90, 129, 116, 112, 110];
    let ordSumLi = seq1.split('').map((char, i) => char.charCodeAt(0) + seq2.charCodeAt(i));
    let loopRegionDict = { "region_pos": [], "region_type": [] };
    let regionTypeFlag = 0;
    let flag = 0;

    for (let i = 0; i < ordSumLi.length; i++) {
        let v = ordSumLi[i];
        if (!matchAscLi.includes(v) && !flag) {
            loopRegionDict["region_pos"].push(i);
            flag = 1;
        }
        if (!gapAscLi.includes(v) && !matchAscLi.includes(v)) {
            regionTypeFlag = 1;
        }
        if (matchAscLi.includes(v) && flag) {
            loopRegionDict["region_pos"].push(i - 1);
            loopRegionDict["region_type"].push(regionTypeFlag);
            flag = 0;
            regionTypeFlag = 0;
        }
    }
    return loopRegionDict;
}

function base2int(base) {
    const trantab = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 };
    return parseInt(base.toUpperCase().split('').map(char => trantab[char]).join(''), 4);
}

function stackEnergy(segment) {
    let deltaH = [-7.9, -8.4, -7.8, -7.2, -8.5, -8.0, -10.6, -7.8, -8.2, -9.8, -8.0, -8.4, -7.2, -8.2, -8.5, -7.9];
    let deltaS = [-22.2, -22.4, -21.0, -20.4, -22.7, -19.9, -27.2, -21.0, -22.2, -24.4, -19.9, -22.4, -21.3, -22.2, -22.7, -22.2];
    let stackDh = 0, stackDs = 0, stackDg = 0;

    for (let i = 0; i < segment.length - 1; i++) {
        let twoMer = segment[i] + segment[i + 1];
        let dIndex = base2int(twoMer);
        stackDh += deltaH[dIndex];
        stackDs += deltaS[dIndex];
    }
    stackDg = stackDh - 310.15 * stackDs / 1000;
    return [stackDh, stackDg];
}

function intermolecularInitiationEnergy() {
    let iiDh = -7.2;
    let iiDg = 1.0;
    return [iiDh, iiDg];
}

function statesCorrectionDg(seq, bulgePos) {
    let stateNum = 1;
    let bulgeBase = seq[bulgePos];
    let forwardIdx = bulgePos, afterwardIdx = bulgePos;

    while (forwardIdx !== 0 && afterwardIdx !== seq.length - 1) {
        if (forwardIdx !== 0 && seq[forwardIdx] === bulgeBase) {
            stateNum += 1;
            forwardIdx += 1;
        } else {
            forwardIdx = 0;
        }

        if (afterwardIdx !== seq.length - 1 && seq[afterwardIdx] === bulgeBase) {
            stateNum += 1;
            afterwardIdx += 1;
        } else {
            afterwardIdx = seq.length - 1;
        }
    }
    return -0.616 * Math.log(stateNum);
}

function bulgeEnergy(segment, bulgeLength, seq, bulgePos) {
    let bulgeLoopDict = {
        "dh": [18.9, -0.6, -2.3, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1, -14.1],
        "dg": [2.9, 2.3, 2.5, 2.7, 3.0, 3.2, 3.4, 3.5, 3.6, 3.7, 3.9, 3.9, 4.0, 4.1, 4.2, 4.3, 4.3, 4.4, 4.4, 4.5, 4.5, 4.6, 4.6, 4.7, 4.7, 4.7, 4.8, 4.9, 4.9, 4.9]
    };

    let bulgeLoopDh = bulgeLength <= 30 ? bulgeLoopDict["dh"][bulgeLength - 1] : bulgeLoopDict["dh"].slice(-1)[0];
    let bulgeLoopDg = bulgeLength <= 30 ? bulgeLoopDict["dg"][bulgeLength - 1] : bulgeLoopDict["dg"].slice(-1)[0];

    if (bulgeLength === 1) {
        let [exStackDh, exStackDg] = stackEnergy(segment[0] + segment[segment.length - 1]);
        let statesDg = statesCorrectionDg(seq, bulgePos);
        bulgeLoopDh += exStackDh;
        bulgeLoopDg += exStackDg + statesDg;
    } else {
        let [atDh, atDg] = ATclosureEnergy(segment);
        bulgeLoopDh += atDh;
        bulgeLoopDg += atDg;

        if (bulgeLength > 6) {
            bulgeLoopDg += 1.75 * 0.616 * Math.log(bulgeLength / 6);
        }
    }
    return [bulgeLoopDh, bulgeLoopDg];
}

function indexIntl11(upstreamBase, downstreamBase, xBase, yBase) {
    let externalIndex = base2int(upstreamBase + downstreamBase);
    let internalIndex = base2int(xBase + yBase);
    return [externalIndex, internalIndex];
}

function indexIntl21(upstreamBase, downstreamBase, xBase, yBase, yNeiborBase) {
    let externalIndex = base2int(yNeiborBase + upstreamBase + downstreamBase);
    let internalIndex = base2int(xBase + yBase);
    return [externalIndex, internalIndex];
}

function indexIntl22(upstreamBase, downstreamBase, xBase1, yBase1, xBase2, yBase2) {
    let externalIndex = base2int(upstreamBase + downstreamBase);
    let internalIndex = base2int(xBase1 + xBase2 + yBase1 + yBase2);
    return [externalIndex, internalIndex];
}

function symmetricIntLoopEnergy(segment1, segment2, loopSum, intLoopEnergyDict) {
    let symmetricIntLoopDh = 0, symmetricIntLoopDg = 0;
    let upstreamBase = segment1[0], downstreamBase = segment1[segment1.length - 1];
    let xBase = segment1[1], yBase = segment2[1];

    if (loopSum === 2) {
        let [externalIndex, internalIndex] = indexIntl11(upstreamBase, downstreamBase, xBase, yBase);
        symmetricIntLoopDh += intLoopEnergyDict["11dh"][externalIndex][internalIndex];
        symmetricIntLoopDg += intLoopEnergyDict["11dg"][externalIndex][internalIndex];
    } else if (loopSum === 3) {
        let yNeiborBase = segment2[2];
        let [externalIndex, internalIndex] = indexIntl21(upstreamBase, downstreamBase, xBase, yBase, yNeiborBase);
        symmetricIntLoopDh += intLoopEnergyDict["21dh"][externalIndex][internalIndex];
        symmetricIntLoopDg += intLoopEnergyDict["21dg"][externalIndex][internalIndex];
    } else {
        let xBase1 = xBase, xBase2 = yBase, yBase1 = segment1[2], yBase2 = segment2[2];
        let [externalIndex, internalIndex] = indexIntl22(upstreamBase, downstreamBase, xBase1, yBase1, xBase2, yBase2);
        symmetricIntLoopDh += intLoopEnergyDict["22dh"][externalIndex][internalIndex];
        symmetricIntLoopDg += intLoopEnergyDict["22dg"][externalIndex][internalIndex];
    }
    console.log("index", externalIndex, internalIndex);
    return [symmetricIntLoopDh, symmetricIntLoopDg];
}


function asymmetryCorrectEnergy(loopDiffAbs) {
    let asymmetryDg = 0.4, asymmetryDh = 0;
    let asymmetryCorrectDh = asymmetryDh * loopDiffAbs;
    let asymmetryCorrectDg = asymmetryDg * loopDiffAbs;
    return [asymmetryCorrectDh, asymmetryCorrectDg];
}

function asymmetricIntLoopInitiationEnergy(loopSum) {
    let initiationDg = [
        3.1, 3.5, 3.9, 4.1, 4.2, 4.3, 4.5, 4.6, 4.6, 4.7,
        4.8, 4.9, 5.0, 5.0, 5.1, 5.1, 5.2, 5.3, 5.3, 5.3,
        5.4, 5.4, 5.5, 5.5, 5.6, 5.6, 5.6
    ];
    let initiationDh = new Array(27).fill(0.0);
    return [initiationDh[loopSum - 4], initiationDg[loopSum - 4]];
}

function asymmetricIntLoopMismatchEnergy(segment1, segment2, intLoopType) {
    let mmDg = [
        // A X, T Y
        [-0.7, -0.3, -0.5, -0.0, -0.6, -0.2, 0.0, -0.3, -0.6, 0.0, -0.4, -0.5, 0.0, -0.3, -0.5, -0.4],
        // C X, G Y
        [-1.0, -0.8, -0.9, 0.0, -0.8, -0.5, 0.0, -0.7, -1.0, 0.0, -0.9, -1.0, 0.0, -0.6, -0.9, -0.9],
        // G X, C Y
        [-1.0, -0.7, -0.8, 0.0, -1.0, -0.6, 0.0, -0.7, -1.0, 0.0, -1.0, -0.8, 0.0, -0.6, -0.9, -0.9],
        // T X, A Y
        [-0.6, -0.4, -0.5, 0.0, -0.5, -0.2, 0.0, -0.5, -0.6, -0.0, -0.4, -0.5, 0.0, -0.3, -0.6, -0.3],
    ];

    let mmDh = [
        // A X, T Y
        [4.0, 5.2, 4.7, 4.3, 4.0, 5.2, 4.7, 4.3, 4.0, 5.2, 4.7, 4.3, 4.0, 5.2, 4.7, 4.3],
        // C X, G Y
        [-4.6, -4.2, -4.7, -5.0, -4.6, -4.2, -4.7, -5.0, -4.6, -4.2, -4.7, -5.0, -4.6, -4.2, -4.7, -5.0],
        // G X, C Y
        [-6.2, -3.3, -5.1, -0.9, -6.2, -3.3, -5.1, -0.9, -6.2, -3.3, -5.1, -0.9, -6.2, -3.3, -5.1, -0.9],
        // T X, A Y
        [2.3, 1.7, 0.7, -5.8, 2.3, 1.7, 0.7, -5.8, 2.3, 1.7, 0.7, -5.8, 2.3, 1.7, 0.7, -5.8],
    ];

    let mismatchDh = 0, mismatchDg = 0;
    let mismatch1Dh = 0, mismatch1Dg = 0;
    let mismatch2Dh = 0, mismatch2Dg = 0;

    if (intLoopType[0] === 1) return [0.0, 0.0];

    let externalIndex = base2int(segment1[0]);
    let internalIndex = base2int(segment1[1] + segment2[1]);
    mismatch1Dh = mmDh[externalIndex][internalIndex];
    mismatch1Dg = mmDg[externalIndex][internalIndex];

    externalIndex = base2int(segment2[segment2.length - 1]);
    internalIndex = base2int(segment2[segment2.length - 2] + segment1[segment1.length - 2]);
    mismatch2Dh = mmDh[externalIndex][internalIndex];
    mismatch2Dg = mmDg[externalIndex][internalIndex];

    mismatchDh = mismatch1Dh + mismatch2Dh;
    mismatchDg = mismatch1Dg + mismatch2Dg;

    return [mismatchDh, mismatchDg];
}

function ATclosureEnergy(segment) {
    let closureAtDh = 0, closureAtDg = 0;
    if (segment[0] === 'A' || segment[0] === 'T') {
        closureAtDh += 3.2;
    }
    if (segment[segment.length - 1] === 'A' || segment[segment.length - 1] === 'T') {
        closureAtDh += 3.2;
    }
    return [closureAtDh, closureAtDg];
}

function asymmetricIntLoopEnergy(segment1, segment2, loopSum, loopDiffAbs, loopType) {
    let asymmetricIntLoopDh = 0, asymmetricIntLoopDg = 0;

    let [liDh, liDg] = asymmetricIntLoopInitiationEnergy(loopSum);
    let [asymmetryCorrectDh, asymmetryCorrectDg] = asymmetryCorrectEnergy(loopDiffAbs);
    let [mmDh, mmDg] = asymmetricIntLoopMismatchEnergy(segment1, segment2, loopType);
    let [atDh, atDg] = ATclosureEnergy(segment1);

    asymmetricIntLoopDh = liDh + asymmetryCorrectDh + mmDh + atDh;
    asymmetricIntLoopDg = liDg + asymmetryCorrectDg + mmDg + atDg;

    return [asymmetricIntLoopDh, asymmetricIntLoopDg];
}

function intLoopEnergy(segment1, segment2, intLoopEnergyDict) {
    let intLoopDh = 0, intLoopDg = 0;

    let firstLoopLength = segment1.length - (segment1.match(/-/g) || []).length - 2;
    let secondLoopLength = segment2.length - (segment2.match(/-/g) || []).length - 2;
    let loopSum = firstLoopLength + secondLoopLength;
    let maxLoopLength = Math.max(firstLoopLength, secondLoopLength);
    let minLoopLength = loopSum - maxLoopLength;
    let loopType = [minLoopLength, maxLoopLength];
    let isSymmetric = loopSum <= 4 && maxLoopLength < 3;

    if (loopSum === 3) {
        if (firstLoopLength === 2) {
            let tmp = segment1.split("").reverse().join("");
            segment1 = segment2.split("").reverse().join("");
            segment2 = tmp;
        }
        if (segment1[segment1.length - 2] !== "-") {
            segment1 = segment1[0] + segment1[segment1.length - 2] + "-" + segment1[segment1.length - 1];
        }
    }

    if (isSymmetric) {
        let [symmetricIntLoopDh, symmetricIntLoopDg] = symmetricIntLoopEnergy(segment1, segment2, loopSum, intLoopEnergyDict);
        intLoopDh += symmetricIntLoopDh;
        intLoopDg += symmetricIntLoopDg;
    } else {
        let loopDiffAbs = Math.abs(firstLoopLength - secondLoopLength);
        let [asymmetricIntLoopDh, asymmetricIntLoopDg] = asymmetricIntLoopEnergy(segment1, segment2, loopSum, loopDiffAbs, loopType);
        intLoopDh += asymmetricIntLoopDh;
        intLoopDg += asymmetricIntLoopDg;
    }

    return [intLoopDh, intLoopDg];
}

function symmetry(seq) {
    let seqLen = seq.length;
    let mp = Math.floor(seqLen / 2);
    if (seqLen % 2 === 1) return false;
    for (let i = 0; i < mp; i++) {
        if (seq[i] !== seq[seqLen - 1 - i]) return false;
    }
    return true;
}

function divalentToMonovalent(divalent, dntp) {
    /*
    Convert divalent salt concentration to monovalent
    */
    if (divalent === 0) dntp = 0;
    if (divalent < dntp) divalent = dntp;
    return 120 * Math.sqrt(divalent - dntp);
}


function calcTmByNN(duplexStr, loopRegionDict) {
    const intLoopEnergyDict = {
        "11dg":[
            // A A
            // T T
            [
                1.3, 2.2, 0.9, 1.0, // AA AC AG AT
                1.7, 2.4, 1.0, 1.4, // CA CC CG CT
                1.7, 2.4, 1.0, 1.4, // GA GC GG GT
                1.7, 2.4, 1.0, 1.4, // TA TC TG TT           
            ],
            // A C
            // T G
            [
                0.8,  1.4, -0.4, 1.0, // AA AC AG AT
                1.6,  2.1,  1.0, 1.6, // CA CC CG CT
                -0.2, 1.0, -1.2, 2.0, // GA GC GG GT
                1.0,  1.4,  2.0, 1.1, // TA TC TG TT    
            ],
            // A G
            // T C
            [
                1.0, 1.7,  0.3, 1.0, // AA AC AG AT
                1.5, 2.0,  1.0, 1.0, // CA CC CG CT
                0.1, 1.0, -0.3, 2.0, // GA GC GG GT
                1.0, 1.4,  2.0, 0.6, // TA TC TG TT 
            ],
            // A T
            // T A
            [
                1.2, 1.7,  0.2, 1.0, // AA AC AG AT
                1.7, 2.7,  1.0, 1.4, // CA CC CG CT
                0.2, 1.0, -0.3, 2.0, // GA GC GG GT
                1.0, 1.4,  2.0, 1.4, // TA TC TG TT 
            ],
            // C A
            // G T
            [
                1.1, 2.1, 0.8, 1.0, // AA AC AG AT
                1.7, 1.8, 1.0, 1.4, // CA CC CG CT
                0.5, 1.0, 0.3, 2.0, // GA GC GG GT
                1.0, 1.4, 2.0, 0.6, // TA TC TG TT 
            ],
            // C C
            // G G
            [
                0.6,  1.2, -0.5, 1.0, // AA AC AG AT
                1.6,  1.5,  1.0, 1.6, // CA CC CG CT
                -0.1, 1.0, -1.3, 2.0, // GA GC GG GT
                1.0,  1.0,  2.0, 0.3, // TA TC TG TT 
            ],
            // C G
            // G C
            [
                0.9, 1.5,  0.1,  1.0, // AA AC AG AT
                1.5, 1.4,  1.0,  1.0, // CA CC CG CT
                0.1, 1.0, -0.3,  2.0, // GA GC GG GT
                1.0, 1.0,  2.0, -0.2, // TA TC TG TT 
            ],
            // C T
            // G A
            [
                1.0, 1.5,  0.1, 1.0, // AA AC AG AT
                1.7, 2.0,  1.0, 1.4, // CA CC CG CT
                0.3, 1.0, -0.3, 2.0, // GA GC GG GT
                1.0, 1.0,  2.0, 0.6, // TA TC TG TT 
            ],
            // G A
            // C T
            [
                0.9,  2.1,  0.5, 1.0, // AA AC AG AT
                1.4,  1.8,  1.0, 1.4, // CA CC CG CT
                -0.1, 1.0, -0.7, 2.0, // GA GC GG GT
                1.0,  2.0,  2.0, 1.1, // TA TC TG TT 
            ],
            // G C
            // C G
            [
                0.3,  1.3, -0.8, 1.0, // AA AC AG AT
                1.3,  1.6,  1.0, 1.6, // CA CC CG CT
                -0.8, 1.0, -2.2, 2.0, // GA GC GG GT
                1.0,  1.6,  2.0, 0.9, // TA TC TG TT 
            ],
            // G G
            // C C
            [
                0.6,  1.6, -0.1, 1.0, // AA AC AG AT
                1.2,  1.5,  1.0, 1.0, // CA CC CG CT
                -0.5, 1.0, -1.3, 2.0, // GA GC GG GT
                1.0,  1.6,  2.0, 0.3, // TA TC TG TT 
            ],
            // G T
            // C A
            [
                0.8,  1.6, -0.2, 1.0, // AA AC AG AT
                1.4,  2.1,  1.0, 1.4, // CA CC CG CT
                -0.4, 1.0, -1.2, 2.0, // GA GC GG GT
                1.0,  1.6,  2.0, 1.1, // TA TC TG TT 
            ],
            // T A
            // A T
            [
                1.4, 2.3, 1.2, 1.0, // AA AC AG AT
                2.3, 2.1, 1.0, 1.7, // CA CC CG CT
                1.2, 1.0, 0.9, 2.0, // GA GC GG GT
                1.0, 1.7, 2.0, 1.4, // TA TC TG TT 
            ],
            // T C
            // A G
            [
                0.9, 1.4, -0.1, 1.0, // AA AC AG AT
                2.1, 1.8,  1.0, 2.0, // CA CC CG CT
                0.5, 1.0, -0.7, 2.0, // GA GC GG GT
                1.0, 1.4,  2.0, 1.1, // TA TC TG TT 
            ],
            // T G
            // A C
            [
                1.1, 1.7, 0.5, 1.0, // AA AC AG AT
                2.1, 1.8, 1.0, 1.4, // CA CC CG CT
                0.8, 1.0, 0.3, 2.0, // GA GC GG GT
                1.0, 1.4, 2.0, 0.6, // TA TC TG TT 
            ],
            // T T
            // A A
            [
                1.3, 1.7, 0.4, 1.0, // AA AC AG AT
                2.2, 2.4, 1.0, 1.7, // CA CC CG CT
                0.9, 1.0, 0.3, 2.0, // GA GC GG GT
                1.0, 1.4, 2.0, 1.4, // TA TC TG TT 
            ]      
    ],
        
        "11dh" :[
            // A A
            // T T
            [
                12.9,1.5,10.7,1.5,
                10.4,14.7,10.4,-1.6,
                1.8,1.5,1.6,1.5,
                10.4,3.5,10.4,-4.8,
            ],
            // A C
            // T G
            [
                10.0,10.9,8.9,10.9,
                4.1,12.8,4.1,1.3,
                4.5,10.9,-5.0,10.9,
                4.1,11.6,4.1,0.5,
            ],
            // A G
            // T C
            [
                5.0,3.9,-6.3,3.9,
                5.3,2.6,5.3,-4.9,
                -2.2,3.9,-9.9,3.9,
                5.3,-12.3,5.3,-7.5,
            ],
            // A T
            // T A
            [
                4.9,9.8,7.1,9.8,
                9.8,7.5,9.8,4.3,
                7.1,9.8,4.1,9.8,
                9.8,4.3,9.8,5.7,
            ],
            // C A
            // G T
            [
                10.0,13.2,11.3,13.2,
                9.6,5.9,9.6,0.1,
                6.9,13.2,-1.1,13.2,
                9.6,3.0,9.6,0.7,
            ],
            // C C
            // G G
            [
                -0.6,4.6,1.3,4.6,
                9.3,7.5,9.3,1.9,
                3.9,4.6,-0.9,4.6,
                9.3,0.5,9.3,2.1,
            ],
            // C G
            // G C
            [
                4.0,5.1,6.1,5.1,
                5.1,4.1,5.1,1.6,
                6.1,5.1,0.1,5.1,
                5.1,1.6,5.1,-2.7,
            ],
            // C T
            // G A
            [
                5.0,5.3,-2.2,5.3,
                3.9,2.6,3.9,-12.3,
                -6.3,5.3,-9.9,5.3,
                3.9,-4.9,3.9,-7.5,
            ],
            // G A
            // C T
            [
                -1.4,13.2,11.3,13.2,
                4.2,12.9,4.2,-0.4,
                6.9,13.2,3.4,13.2,
                4.2,3.7,4.2,-10.0,
            ],
            // G C
            // C G
            [
                2.6,3.9,2.4,3.9,
                3.9,6.0,3.9,-0.6,
                2.4,3.9,-2.8,3.9,
                3.9,-0.6,3.9,0.3,
            ],
            // G G
            // C C
            [
                -0.6,9.3,3.9,9.3,
                4.6,7.5,4.6,0.5,
                1.3,9.3,-0.9,9.3,
                4.6,1.9,4.6,2.1,
            ],
            // G T
            // C A
            [
                10.0,4.1,4.5,4.1,
                10.9,12.8,10.9,11.6,
                8.9,4.1,-5.0,4.1,
                10.9,1.3,10.9,0.5,
            ],
            // T A
            // A T
            [
                12.1,14.3,11.9,14.3,
                14.3,17.4,14.3,2.5,
                11.9,14.3,8.4,14.3,
                14.3,2.5,14.3,6.2,
            ],
            // T C
            // A G
            [
                -1.4,4.2,6.9,4.2,
                13.2,12.9,13.2,3.7,
                11.3,4.2,3.4,4.2,
                13.2,-0.4,13.2,-10.0,
            ],
            // T G
            // A C
            [
                10.0,9.6,6.9,9.6,
                13.2,5.9,13.2,3.0,
                11.3,9.6,-1.1,9.6,
                13.2,0.1,13.2,0.7,
            ],
            // T T
            // A A
            [
                12.9,10.4,1.8,10.4,
                1.5,14.7,1.5,3.5,
                10.7,10.4,1.6,10.4,
                1.5,-1.6,1.5,-4.8,
            ],
    ],
        "21dg" : [
            // A  A
            // T AT
            [
                1.8,2.0,1.9,2.5,
                2.1,2.3,2.2,2.2,
                1.9,2.3,1.9,2.5,
                3.1,2.7,3.1,2.5,
            ],
            // A  C
            // T AG
            [
                1.3,1.5,1.4,2.0,
                1.8,2.0,1.9,1.9,
                1.5,1.9,1.5,2.1,
                2.6,2.2,2.6,2.0,
            ],
            // A  G
            // T AC
            [
                1.4,1.6,1.5,2.1,
                1.9,2.1,2.0,2.0,
                1.6,2.0,1.6,2.2,
                2.7,2.3,2.7,2.1,
            ],
            // A  T
            // T AA
            [
                1.7,1.9,1.8,2.4,
                2.1,2.3,2.2,2.2,
                1.9,2.3,1.9,2.5,
                3.1,2.7,3.1,2.5,
            ],
            // A  A
            // T CT
            [
                1.9,2.1,2.0,2.6,
                2.3,2.5,2.4,2.4,
                2.0,2.4,2.0,2.6,
                2.7,2.3,2.7,2.1,
            ],
            // A  C
            // T CG
            [
                1.6,1.8,1.7,2.3,
                2.0,2.2,2.1,2.1,
                1.2,1.6,1.2,1.8,
                2.4,2.0,2.4,1.8,
            ],
            // A  G
            // T CC
            [
                1.7,1.9,1.8,2.4,
                2.1,2.3,2.2,2.2,
                1.7,2.1,1.7,2.3,
                2.4,2.0,2.4,1.8,
            ],
            // A  T
            // T CA
            [
                1.9,2.1,2.0,2.6,
                2.3,2.5,2.4,2.4,
                2.1,2.5,2.1,2.7,
                2.7,2.3,2.7,2.1,
            ],
            // A  A
            // T GT
            [
                1.8,2.0,1.9,2.5,
                2.6,2.8,2.7,2.7,
                1.9,2.3,1.9,2.5,
                3.1,2.7,3.1,2.5,
            ],
            // A  C
            // T GG
            [
                1.4,1.6,1.5,2.1,
                1.8,2.0,1.9,1.9,
                1.5,1.9,1.5,2.1,
                2.8,2.4,2.8,2.2,
            ],
            // A  G
            // T GC
            [
                1.5,1.7,1.6,2.2,
                1.5,1.7,1.6,1.6,
                1.6,2.0,1.6,2.2,
                2.6,2.2,2.6,2.0,
            ],
            // A  T
            // T GA
            [
                1.8,2.0,1.9,2.5,
                2.4,2.6,2.5,2.5,
                1.9,2.3,1.9,2.5,
                3.1,2.7,3.1,2.5,
            ],
            // A  A
            // T TT
            [
                2.4,2.6,2.5,3.1,
                2.2,2.4,2.3,2.3,
                2.4,2.8,2.4,3.0,
                2.5,2.1,2.5,1.9,
            ],
            // A  C
            // T TG
            [
                2.0,2.2,2.1,2.7,
                1.9,2.1,2.0,2.0,
                2.1,2.5,2.1,2.7,
                2.2,1.8,2.2,1.6,
            ],
            // A  G
            // T TC
            [
                2.2,2.4,2.3,2.9,
                1.9,2.1,2.0,2.0,
                2.1,2.5,2.1,2.7,
                2.2,1.8,2.2,1.6,
            ],
            // A  T
            // T TA
            [
                2.4,2.6,2.5,3.1,
                2.2,2.4,2.3,2.3,
                2.5,2.9,2.5,3.1,
                2.5,2.1,2.5,1.9,
            ],
            // C  A
            // G AT
            [
                1.5,1.8,1.6,2.1,
                1.9,2.1,1.8,1.9,
                1.6,1.4,1.6,2.0,
                2.9,2.4,2.7,2.2,
            ],
            // C  C
            // G AG
            [
                1.0,1.3,1.1,1.6,
                1.6,1.8,1.5,1.6,
                1.2,1.0,1.2,1.6,
                2.4,1.9,2.2,1.7,
            ],
            // C  G
            // G AC
            [
                1.1,1.4,1.2,1.7,
                1.7,1.9,1.6,1.7,
                1.3,1.1,1.3,1.7,
                2.5,2.0,2.3,1.8,
            ],
            // C  T
            // G AA
            [
                1.4,1.7,1.5,2.0,
                1.9,2.1,1.8,1.9,
                1.6,1.4,1.6,2.0,
                2.9,2.4,2.7,2.2,
            ],
            // C  A
            // G CT
            [
                1.6,1.9,1.7,2.2,
                2.1,2.3,2.0,2.1,
                1.7,1.5,1.7,2.1,
                2.5,2.0,2.3,1.8,
            ],
            // C  C
            // G CG
            [
                1.3,1.6,1.4,1.9,
                1.8,2.0,1.7,1.8,
                0.9,0.7,0.9,1.3,
                2.2,1.7,2.0,1.5,
            ],
            // C  G
            // G CC
            [
                1.4,1.7,1.5,2.0,
                1.9,2.1,1.8,1.9,
                1.4,1.2,1.4,1.8,
                2.2,1.7,2.0,1.5,
            ],
            // C  T
            // G CA
            [
                1.6,1.9,1.7,2.2,
                2.1,2.3,2.0,2.1,
                1.8,1.6,1.8,2.2,
                2.5,2.0,2.3,1.8,
            ],
            // C  A
            // G GT
            [
                1.5,1.8,1.6,2.1,
                2.4,2.6,2.3,2.4,
                1.6,1.4,1.6,2.0,
                2.9,2.4,2.7,2.2,
            ],
            // C  C
            // G GG
            [
                1.1,1.4,1.2,1.7,
                1.6,1.8,1.5,1.6,
                1.2,1.0,1.2,1.6,
                2.6,2.1,2.4,1.9,
            ],
            // C  G
            // G GC
            [
                1.2,1.5,1.3,1.8,
                1.3,1.5,1.2,1.3,
                1.3,1.1,1.3,1.7,
                2.4,1.9,2.2,1.7,
            ],
            // C  T
            // G GA
            [
                1.5,1.8,1.6,2.1,
                2.2,2.4,2.1,2.2,
                1.6,1.4,1.6,2.0,
                2.9,2.4,2.7,2.2,
            ],
            // C  A
            // G TT
            [
                2.1,2.4,2.2,2.7,
                2.0,2.2,1.9,2.0,
                2.1,1.9,2.1,2.5,
                2.3,1.8,2.1,1.6,
            ],
            // C  C
            // G TG
            [
                1.7,2.0,1.8,2.3,
                1.7,1.9,1.6,1.7,
                1.8,1.6,1.8,2.2,
                2.0,1.5,1.8,1.3,
            ],
            // C  G
            // G TC
            [
                1.9,2.2,2.0,2.5,
                1.7,1.9,1.6,1.7,
                1.8,1.6,1.8,2.2,
                2.0,1.5,1.8,1.3,
            ],
            // C  T
            // G TA
            [
                2.1,2.4,2.2,2.7,
                2.0,2.2,1.9,2.0,
                2.2,2.0,2.2,2.6,
                2.3,1.8,2.1,1.6,
            ],
            // G  A
            // C AT
            [
                1.4,1.7,1.5,2.0,
                1.8,2.0,1.3,1.9,
                1.5,1.7,1.5,2.2,
                2.7,2.4,2.7,2.2,
            ],
            // G  C
            // C AG
            [
                0.9,1.2,1.0,1.5,
                1.5,1.7,1.0,1.6,
                1.1,1.3,1.1,1.8,
                2.2,1.9,2.2,1.7,
            ],
            // G  G
            // C AC
            [
                1.0,1.3,1.1,1.6,
                1.6,1.8,1.1,1.7,
                1.2,1.4,1.2,1.9,
                2.3,2.0,2.3,1.8,
            ],
            // G  T
            // C AA
            [
                1.3,1.6,1.4,1.9,
                1.8,2.0,1.3,1.9,
                1.5,1.7,1.5,2.2,
                2.7,2.4,2.7,2.2,
            ],
            // G  A
            // C CT
            [
                1.5,1.8,1.6,2.1,
                2.0,2.2,1.5,2.1,
                1.6,1.8,1.6,2.3,
                2.3,2.0,2.3,1.8,
            ],
            // G  C
            // C CG
            [
                1.2,1.5,1.3,1.8,
                1.7,1.9,1.2,1.8,
                0.8,1.0,0.8,1.5,
                2.0,1.7,2.0,1.5,
            ],
            // G  G
            // C CC
            [
                1.3,1.6,1.4,1.9,
                1.8,2.0,1.3,1.9,
                1.3,1.5,1.3,2.0,
                2.0,1.7,2.0,1.5,
            ],
            // G  T
            // C CA
            [
                1.5,1.8,1.6,2.1,
                2.0,2.2,1.5,2.1,
                1.7,1.9,1.7,2.4,
                2.3,2.0,2.3,1.8,
            ],
            // G  A
            // C GT
            [
                1.4,1.7,1.5,2.0,
                2.3,2.5,1.8,2.4,
                1.5,1.7,1.5,2.2,
                2.7,2.4,2.7,2.2,
            ],
            // G  C
            // C GG
            [
                1.0,1.3,1.1,1.6,
                1.5,1.7,1.0,1.6,
                1.1,1.3,1.1,1.8,
                2.4,2.1,2.4,1.9,
            ],
            // G  G
            // C GC
            [
                1.1,1.4,1.2,1.7,
                1.2,1.4,0.7,1.3,
                1.2,1.4,1.2,1.9,
                2.2,1.9,2.2,1.7,
            ],
            // G  T
            // C GA
            [
                1.4,1.7,1.5,2.0,
                2.1,2.3,1.6,2.2,
                1.5,1.7,1.5,2.2,
                2.7,2.4,2.7,2.2,
            ],
            // G  A
            // C TT
            [
                2.0,2.3,2.1,2.6,
                1.9,2.1,1.4,2.0,
                2.0,2.2,2.0,2.7,
                2.1,1.8,2.1,1.6,
            ],
            // G  C
            // C TG
            [
                1.6,1.9,1.7,2.2,
                1.6,1.8,1.1,1.7,
                1.7,1.9,1.7,2.4,
                1.8,1.5,1.8,1.3,
            ],
            // G  G
            // C TC
            [
                1.8,2.1,1.9,2.4,
                1.6,1.8,1.1,1.7,
                1.7,1.9,1.7,2.4,
                1.8,1.5,1.8,1.3,
            ],
            // G  T
            // C TA
            [
                2.0,2.3,2.1,2.6,
                1.9,2.1,1.4,2.0,
                2.1,2.3,2.1,2.8,
                2.1,1.8,2.1,1.6,
            ],
            // T  A
            // A AT
            [
                1.9,2.0,1.9,2.5,
                2.1,2.3,2.1,2.2,
                1.9,2.5,1.9,2.5,
                3.1,2.7,3.0,2.5,
            ],
            // T  C
            // A AG
            [
                1.4,1.5,1.4,2.0,
                1.8,2.0,1.8,1.9,
                1.5,2.1,1.5,2.1,
                2.6,2.2,2.5,2.0,
            ],
            // T  G
            // A AC
            [
                1.5,1.6,1.5,2.1,
                1.9,2.1,1.9,2.0,
                1.6,2.2,1.6,2.2,
                2.7,2.3,2.6,2.1,
            ],
            // T  T
            // A AA
            [
                1.8,1.9,1.8,2.4,
                2.1,2.3,2.1,2.2,
                1.9,2.5,1.9,2.5,
                3.1,2.7,3.0,2.5,
            ],
            // T  A
            // A CT
            [
                2.0,2.1,2.0,2.6,
                2.3,2.5,2.3,2.4,
                2.0,2.6,2.0,2.6,
                2.7,2.3,2.6,2.1,
            ],
            // T  C
            // A CG
            [
                1.7,1.8,1.7,2.3,
                2.0,2.2,2.0,2.1,
                1.2,1.8,1.2,1.8,
                2.4,2.0,2.3,1.8,
            ],
            // T  G
            // A CC
            [
                1.8,1.9,1.8,2.4,
                2.1,2.3,2.1,2.2,
                1.7,2.3,1.7,2.3,
                2.4,2.0,2.3,1.8,
            ],
            // T  T
            // A CA
            [
                2.0,2.1,2.0,2.6,
                2.3,2.5,2.3,2.4,
                2.1,2.7,2.1,2.7,
                2.7,2.3,2.6,2.1,
            ],
            // T  A
            // A GT
            [
                1.9,2.0,1.9,2.5,
                2.6,2.8,2.6,2.7,
                1.9,2.5,1.9,2.5,
                3.1,2.7,3.0,2.5,
            ],
            // T  C
            // A GG
            [
                1.5,1.6,1.5,2.1,
                1.8,2.0,1.8,1.9,
                1.5,2.1,1.5,2.1,
                2.8,2.4,2.7,2.2,
            ],
            // T  G
            // A GC
            [
                1.6,1.7,1.6,2.2,
                1.5,1.7,1.5,1.6,
                1.6,2.2,1.6,2.2,
                2.6,2.2,2.5,2.0,
            ],
            // T  T
            // A GA
            [
                1.9,2.0,1.9,2.5,
                2.4,2.6,2.4,2.5,
                1.9,2.5,1.9,2.5,
                3.1,2.7,3.0,2.5,
            ],
            // T  A
            // A TT
            [
                2.5,2.6,2.5,3.1,
                2.2,2.4,2.2,2.3,
                2.4,3.0,2.4,3.0,
                2.5,2.1,2.4,1.9,
            ],
            // T  C
            // A TG
            [
                2.1,2.2,2.1,2.7,
                1.9,2.1,1.9,2.0,
                2.1,2.7,2.1,2.7,
                2.2,1.8,2.1,1.6,
            ],
            // T  G
            // A TC
            [
                2.3,2.4,2.3,2.9,
                1.9,2.1,1.9,2.0,
                2.1,2.7,2.1,2.7,
                2.2,1.8,2.1,1.6,
            ],
            // T  T
            // A TA
            [
                2.5,2.6,2.5,3.1,
                2.2,2.4,2.2,2.3,
                2.5,3.1,2.5,3.1,
                2.5,2.1,2.4,1.9,
            ],
    ],
        "21dh" : [
            // A  A
            // T AT
            [
                31.8,31.8,29.6,31.8,
                29.3,29.3,29.3,29.3,
                20.7,20.7,20.7,20.7,
                29.3,29.3,29.3,14.1,
            ],
            // A  C
            // T AG
            [
                28.9,29.8,27.8,29.8,
                23.0,23.0,23.0,20.2,
                23.4,23.4,13.9,23.4,
                23.0,23.0,23.0,19.4,
            ],
            // A  G
            // T AC
            [
                23.9,23.9,12.6,23.9,
                24.2,24.2,24.2,14.0,
                16.7,16.7,16.7,16.7,
                24.2,6.6,24.2,11.4,
            ],
            // A  T
            // T AA
            [
                23.8,28.7,26.0,23.8,
                28.7,28.7,28.7,23.2,
                26.0,26.0,26.0,26.0,
                28.7,23.2,28.7,28.7,
            ],
            // A  A
            // T CT
            [
                31.8,20.4,29.6,20.4,
                29.3,33.6,29.3,17.3,
                20.7,20.4,20.5,20.4,
                29.3,22.4,29.3,14.1,
            ],
            // A  C
            // T CG
            [
                29.8,29.8,27.8,29.8,
                23.0,31.7,23.0,20.2,
                23.4,29.8,13.9,29.8,
                23.0,30.5,23.0,19.4,
            ],
            // A  G
            // T CC
            [
                22.8,22.8,12.6,22.8,
                24.2,21.5,24.2,14.0,
                16.7,22.8,9.0,22.8,
                6.6,6.6,6.6,6.6,
            ],
            // A  T
            // T CA
            [
                28.7,28.7,26.0,28.7,
                28.7,26.4,28.7,23.2,
                26.0,28.7,23.0,28.7,
                23.2,23.2,23.2,23.2,
            ],
            // A  A
            // T GT
            [
                29.6,29.6,29.6,29.6,
                29.3,29.3,29.3,29.3,
                20.7,20.5,20.5,20.5,
                29.3,29.3,29.3,14.1,
            ],
            // A  C
            // T GG
            [
                27.8,27.8,27.8,27.8,
                23.0,23.0,23.0,20.2,
                13.9,13.9,13.9,13.9,
                23.0,23.0,23.0,19.4,
            ],
            // A  G
            // T GC
            [
                12.6,12.6,12.6,12.6,
                24.2,24.2,24.2,14.0,
                16.7,9.0,9.0,9.0,
                24.2,6.6,24.2,11.4,
            ],
            // A  T
            // T GA
            [
                26.0,26.0,26.0,26.0,
                28.7,28.7,28.7,23.2,
                26.0,23.0,23.0,23.0,
                28.7,23.2,28.7,28.7,
            ],
            // A  A
            // T TT
            [
                31.8,20.4,29.6,20.4,
                17.3,17.3,17.3,17.3,
                20.7,20.4,20.5,20.4,
                14.1,14.1,14.1,14.1,
            ],
            // A  C
            // T TG
            [
                29.8,29.8,27.8,29.8,
                20.2,20.2,20.2,20.2,
                23.4,29.8,13.9,29.8,
                19.4,19.4,19.4,19.4,
            ],
            // A  G
            // T TC
            [
                22.8,22.8,12.6,22.8,
                14.0,14.0,14.0,14.0,
                16.7,22.8,9.0,22.8,
                11.4,6.6,11.4,11.4,
            ],
            // A  T
            // T TA
            [
                28.7,28.7,26.0,28.7,
                23.2,23.2,23.2,23.2,
                26.0,28.7,23.0,28.7,
                28.7,23.2,28.7,24.6,
            ],
            // C  A
            // G AT
            [
                28.9,28.9,30.2,28.9,
                28.5,28.5,28.5,19.0,
                25.8,25.8,17.8,25.8,
                28.5,21.9,28.5,19.6,
            ],
            // C  C
            // G AG
            [
                18.3,23.5,20.2,23.5,
                28.2,28.2,28.2,20.8,
                22.8,22.8,18.0,22.8,
                28.2,19.4,28.2,21.0,
            ],
            // C  G
            // G AC
            [
                22.9,22.9,25.0,22.9,
                24.0,24.0,24.0,20.5,
                25.0,25.0,25.0,25.0,
                24.0,20.5,24.0,16.2,
            ],
            // C  T
            // G AA
            [
                23.9,24.2,16.7,24.2,
                22.8,22.8,22.8,6.6,
                12.6,12.6,9.0,12.6,
                22.8,14.0,22.8,11.4,
            ],
            // C  A
            // G CT
            [
                28.9,32.1,30.2,32.1,
                28.5,24.8,28.5,19.0,
                25.8,32.1,17.8,32.1,
                21.9,21.9,21.9,19.6,
            ],
            // C  C
            // G CG
            [
                23.5,23.5,20.2,23.5,
                28.2,26.4,28.2,20.8,
                22.8,23.5,18.0,23.5,
                19.4,19.4,19.4,19.4,
            ],
            // C  G
            // G CC
            [
                22.9,24.0,25.0,24.0,
                24.0,23.0,24.0,20.5,
                25.0,24.0,19.0,24.0,
                20.5,20.5,20.5,16.2,
            ],
            // C  T
            // G CA
            [
                24.2,24.2,16.7,24.2,
                22.8,21.5,22.8,6.6,
                12.6,24.2,9.0,24.2,
                14.0,14.0,14.0,11.4,
            ],
            // C  A
            // G GT
            [
                30.2,30.2,30.2,30.2,
                28.5,28.5,28.5,19.0,
                17.8,17.8,17.8,17.8,
                28.5,21.9,28.5,19.6,
            ],
            // C  C
            // G GG
            [
                20.2,20.2,20.2,20.2,
                28.2,28.2,28.2,20.8,
                18.0,18.0,18.0,18.0,
                28.2,19.4,28.2,21.0,
            ],
            // C  G
            // G GC
            [
                25.0,25.0,25.0,25.0,
                24.0,24.0,24.0,20.5,
                25.0,19.0,19.0,19.0,
                24.0,20.5,24.0,16.2,
            ],
            // C  T
            // G GA
            [
                16.7,16.7,16.7,16.7,
                22.8,22.8,22.8,6.6,
                9.0,9.0,9.0,9.0,
                22.8,14.0,22.8,11.4,
            ],
            // C  A
            // G TT
            [
                28.9,32.1,30.2,32.1,
                19.0,19.0,19.0,19.0,
                25.8,32.1,17.8,32.1,
                19.6,19.6,19.6,19.6,
            ],
            // C  C
            // G TG
            [
                23.5,23.5,20.2,23.5,
                20.8,20.8,20.8,20.8,
                22.8,23.5,18.0,23.5,
                21.0,21.0,21.0,21.0,
            ],
            // C  G
            // G TC
            [
                22.9,24.0,25.0,24.0,
                20.5,20.5,20.5,20.5,
                25.0,24.0,19.0,24.0,
                16.2,16.2,16.2,16.2,
            ],
            // C  T
            // G TA
            [
                24.2,24.2,16.7,24.2,
                6.6,6.6,6.6,6.6,
                12.6,24.2,9.0,24.2,
                11.4,11.4,11.4,11.4,
            ],
            // G  A
            // C AT
            [
                17.5,17.5,17.5,17.5,
                23.1,23.1,23.1,23.1,
                25.8,25.8,25.8,25.8,
                23.1,23.1,23.1,8.9,
            ],
            // G  C
            // C AG
            [
                21.5,21.5,21.3,21.5,
                22.8,22.8,22.8,22.8,
                21.3,21.3,16.1,21.3,
                22.8,22.8,22.8,22.8,
            ],
            // G  G
            // C AC
            [
                18.3,18.3,22.8,18.3,
                23.5,23.5,23.5,23.5,
                20.2,20.2,20.2,20.2,
                23.5,23.5,23.5,23.5,
            ],
            // G  T
            // C AA
            [
                28.9,23.0,23.4,23.0,
                29.8,29.8,29.8,29.8,
                27.8,27.8,13.9,27.8,
                29.8,20.2,29.8,19.4,
            ],
            // G  A
            // C CT
            [
                17.5,32.1,30.2,32.1,
                23.1,31.8,23.1,18.5,
                25.8,32.1,22.3,32.1,
                23.1,22.6,23.1,8.9,
            ],
            // G  C
            // C CG
            [
                21.5,22.8,21.3,22.8,
                22.8,24.9,22.8,18.3,
                21.3,22.8,16.1,22.8,
                22.8,18.3,22.8,18.3,
            ],
            // G  G
            // C CC
            [
                18.3,28.2,22.8,28.2,
                23.5,26.4,23.5,19.4,
                20.2,28.2,18.0,28.2,
                23.5,20.8,23.5,21.0,
            ],
            // G  T
            // C CA
            [
                23.0,23.0,23.4,23.0,
                29.8,31.7,29.8,30.5,
                27.8,23.0,13.9,23.0,
                20.2,20.2,20.2,20.2,
            ],
            // G  A
            // C GT
            [
                30.2,30.2,30.2,30.2,
                23.1,23.1,23.1,23.1,
                25.8,22.3,22.3,22.3,
                23.1,23.1,23.1,8.9,
            ],
            // G  C
            // C GG
            [
                21.3,21.3,21.3,21.3,
                22.8,22.8,22.8,22.8,
                16.1,16.1,16.1,16.1,
                22.8,22.8,22.8,22.8,
            ],
            // G  G
            // C GC
            [
                22.8,22.8,22.8,22.8,
                23.5,23.5,23.5,23.5,
                20.2,18.0,18.0,18.0,
                23.5,23.5,23.5,23.5,
            ],
            // G  T
            // C GA
            [
                23.4,23.4,23.4,23.4,
                29.8,29.8,29.8,29.8,
                13.9,13.9,13.9,13.9,
                29.8,20.2,29.8,19.4,
            ],
            // G  A
            // C TT
            [
                17.5,32.1,30.2,32.1,
                23.1,18.5,23.1,18.5,
                25.8,32.1,22.3,32.1,
                8.9,8.9,8.9,8.9,
            ],
            // G  C
            // C TG
            [
                21.5,22.8,21.3,22.8,
                22.8,18.3,22.8,18.3,
                21.3,22.8,16.1,22.8,
                22.8,18.3,22.8,19.2,
            ],
            // G  G
            // C TC
            [
                18.3,28.2,22.8,28.2,
                23.5,19.4,23.5,19.4,
                20.2,28.2,18.0,28.2,
                23.5,21.0,23.5,21.0,
            ],
            // G  T
            // C TA
            [
                23.0,23.0,23.4,23.0,
                29.8,30.5,29.8,30.5,
                27.8,23.0,13.9,23.0,
                19.4,19.4,19.4,19.4,
            ],
            // T  A
            // A AT
            [
                31.0,31.0,30.8,31.0,
                33.2,33.2,33.2,21.4,
                30.8,30.8,30.8,30.8,
                33.2,21.4,33.2,25.1,
            ],
            // T  C
            // A AG
            [
                17.5,17.5,25.8,17.5,
                32.1,32.1,32.1,32.1,
                30.2,30.2,22.3,30.2,
                32.1,18.5,32.1,8.9,
            ],
            // T  G
            // A AC
            [
                28.9,28.9,25.8,28.9,
                32.1,32.1,32.1,21.9,
                30.2,30.2,17.8,30.2,
                32.1,19.0,32.1,19.6,
            ],
            // T  T
            // A AA
            [
                31.8,29.3,20.7,29.3,
                20.4,20.4,20.4,22.4,
                29.6,29.6,20.5,29.6,
                20.4,17.3,20.4,14.1,
            ],
            // T  A
            // A CT
            [
                31.0,33.2,30.8,33.2,
                33.2,36.3,33.2,21.4,
                30.8,33.2,27.3,33.2,
                21.4,21.4,21.4,21.4,
            ],
            // T  C
            // A CG
            [
                17.5,23.1,25.8,23.1,
                32.1,31.8,32.1,31.8,
                30.2,23.1,22.3,23.1,
                18.5,18.5,18.5,8.9,
            ],
            // T  G
            // A CC
            [
                28.5,28.5,25.8,28.5,
                32.1,24.8,32.1,21.9,
                30.2,28.5,17.8,28.5,
                19.0,19.0,19.0,19.6,
            ],
            // T  T
            // A CA
            [
                29.3,29.3,20.7,29.3,
                20.4,33.6,20.4,22.4,
                29.6,29.3,20.5,29.3,
                17.3,17.3,17.3,14.1,
            ],
            // T  A
            // A GT
            [
                30.8,30.8,30.8,30.8,
                33.2,33.2,33.2,21.4,
                30.8,27.3,27.3,27.3,
                33.2,21.4,33.2,25.1,
            ],
            // T  C
            // A GG
            [
                25.8,25.8,25.8,25.8,
                32.1,32.1,32.1,32.1,
                22.3,22.3,22.3,22.3,
                32.1,18.5,32.1,8.9,
            ],
            // T  G
            // A GC
            [
                25.8,25.8,25.8,25.8,
                32.1,32.1,32.1,21.9,
                17.8,17.8,17.8,17.8,
                32.1,19.0,32.1,19.6,
            ],
            // T  T
            // A GA
            [
                20.7,20.7,20.7,20.7,
                20.4,20.4,20.4,22.4,
                20.5,20.5,20.5,20.5,
                20.4,17.3,20.4,14.1,
            ],
            // T  A
            // A TT
            [
                31.0,33.2,30.8,33.2,
                21.4,21.4,21.4,21.4,
                30.8,33.2,27.3,33.2,
                25.1,21.4,25.1,25.1,
            ],
            // T  C
            // A TG
            [
                17.5,23.1,25.8,23.1,
                32.1,22.6,32.1,22.6,
                30.2,23.1,22.3,23.1,
                8.9,8.9,8.9,8.9,
            ],
            // T  G
            // A TC
            [
                28.5,28.5,25.8,28.5,
                21.9,21.9,21.9,21.9,
                30.2,28.5,17.8,28.5,
                19.6,19.6,19.6,19.6,
            ],
            // T  T
            // A TA
            [
                29.3,29.3,20.7,29.3,
                22.4,22.4,22.4,22.4,
                29.6,29.3,20.5,29.3,
                14.1,14.1,14.1,14.1,
            ],
    ],
        "22dg" : [
            // A A
            // T T
            [
                1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
                2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
                2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
                2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
                2.2,2.3,2.2,2.8,2.3,2.5,2.8,2.4,2.2,2.3,2.2,2.7,2.8,2.4,2.8,2.2,
                2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
                2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
                2.3,2.4,2.3,2.9,2.4,2.6,2.9,2.5,2.3,2.4,2.3,2.8,2.9,2.5,2.9,2.3,
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
                2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
                2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
                2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
                2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
            ],
            // A C
            // T G
            [
                1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
                1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
                2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
                1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
                1.7,2.0,1.8,2.4,2.0,2.2,2.0,2.1,1.8,1.5,1.8,2.4,2.3,2.1,2.5,1.9,
                1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
                1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
                1.8,2.1,1.9,2.5,2.1,2.3,2.1,2.2,1.9,1.6,1.9,2.5,2.4,2.2,2.6,2.0,
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
                2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
                2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
                1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
                2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
            ],
            // A G
            // T C
            [
                1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
                1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
                2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
                1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
                1.8,2.1,1.9,2.6,2.1,2.3,1.7,2.1,1.9,2.0,1.9,2.4,2.4,2.1,2.3,1.9,
                1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
                1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
                1.9,2.2,2.0,2.7,2.2,2.4,1.8,2.2,2.0,2.1,2.0,2.5,2.5,2.2,2.4,2.0,
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
                2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
                2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
                1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
                2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
            ],
            // A T
            // T A
            [
                1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
                1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
                2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
                1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
                2.1,2.3,2.2,2.8,2.3,2.5,2.6,2.4,2.2,2.4,2.2,2.8,2.8,2.4,2.8,2.2,
                2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
                2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
                2.2,2.4,2.3,2.9,2.4,2.6,2.7,2.5,2.3,2.5,2.3,2.9,2.9,2.5,2.9,2.3,
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
                2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
                2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
                2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
                2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
            ],
            // C A
            // G T
            [
                1.5,1.6,1.5,2.1,1.6,1.8,2.1,1.7,1.5,1.6,1.5,2.0,2.1,1.7,2.1,1.5,
                1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
                1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
                2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
                1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
                2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
                1.7,1.8,1.7,2.3,1.8,2.0,2.3,1.9,1.7,1.8,1.7,2.2,2.3,1.9,2.3,1.7,
                1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
                1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
                1.4,1.5,1.4,2.0,1.5,1.7,2.0,1.6,1.4,1.5,1.4,1.9,2.0,1.6,2.0,1.4,
                1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
                2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
                2.3,2.4,2.3,2.9,2.4,2.6,2.9,2.5,2.3,2.4,2.3,2.8,2.9,2.5,2.9,2.3,
                1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
                2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
                1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
            ],
            // C C
            // G G
            [
                1.0,1.3,1.1,1.7,1.3,1.5,1.3,1.4,1.1,0.8,1.1,1.7,1.6,1.4,1.8,1.2,
                1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
                1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
                1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
                1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
                1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
                1.2,1.5,1.3,1.9,1.5,1.7,1.5,1.6,1.3,1.0,1.3,1.9,1.8,1.6,2.0,1.4,
                1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
                1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
                0.9,1.2,1.0,1.6,1.2,1.4,1.2,1.3,1.0,0.7,1.0,1.6,1.5,1.3,1.7,1.1,
                1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
                1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
                1.8,2.1,1.9,2.5,2.1,2.3,2.1,2.2,1.9,1.6,1.9,2.5,2.4,2.2,2.6,2.0,
                1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
                1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
                1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
            ],
            // C G
            // G C
            [
                1.1,1.4,1.2,1.9,1.4,1.6,1.0,1.4,1.2,1.3,1.2,1.7,1.7,1.4,1.6,1.2,
                1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
                1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
                1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
                1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
                1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
                1.3,1.6,1.4,2.1,1.6,1.8,1.2,1.6,1.4,1.5,1.4,1.9,1.9,1.6,1.8,1.4,
                1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
                1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
                1.0,1.3,1.1,1.8,1.3,1.5,0.9,1.3,1.1,1.2,1.1,1.6,1.6,1.3,1.5,1.1,
                1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
                1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
                1.9,2.2,2.0,2.7,2.2,2.4,1.8,2.2,2.0,2.1,2.0,2.5,2.5,2.2,2.4,2.0,
                1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
                1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
                1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
            ],
            // C T
            // G A
            [
                1.4,1.6,1.5,2.1,1.6,1.8,1.9,1.7,1.5,1.7,1.5,2.1,2.1,1.7,2.1,1.5,
                1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
                1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
                2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
                1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
                1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
                1.6,1.8,1.7,2.3,1.8,2.0,2.1,1.9,1.7,1.9,1.7,2.3,2.3,1.9,2.3,1.7,
                1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
                1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
                1.3,1.5,1.4,2.0,1.5,1.7,1.8,1.6,1.4,1.6,1.4,2.0,2.0,1.6,2.0,1.4,
                1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
                1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
                2.2,2.4,2.3,2.9,2.4,2.6,2.7,2.5,2.3,2.5,2.3,2.9,2.9,2.5,2.9,2.3,
                1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
                2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
                1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
            ],
            // G A
            // C T
            [
                1.4,1.5,1.4,2.0,1.5,1.7,2.0,1.6,1.4,1.5,1.4,1.9,2.0,1.6,2.0,1.4,
                1.7,1.8,1.7,2.3,1.8,2.0,2.3,1.9,1.7,1.8,1.7,2.2,2.3,1.9,2.3,1.7,
                1.5,1.6,1.5,2.1,1.6,1.8,2.1,1.7,1.5,1.6,1.5,2.0,2.1,1.7,2.1,1.5,
                2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
                1.7,1.8,1.7,2.3,1.8,2.0,2.3,1.9,1.7,1.8,1.7,2.2,2.3,1.9,2.3,1.7,
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
                1.2,1.3,1.2,1.8,1.3,1.5,1.8,1.4,1.2,1.3,1.2,1.7,1.8,1.4,1.8,1.2,
                1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
                1.5,1.6,1.5,2.1,1.6,1.8,2.1,1.7,1.5,1.6,1.5,2.0,2.1,1.7,2.1,1.5,
                1.7,1.8,1.7,2.3,1.8,2.0,2.3,1.9,1.7,1.8,1.7,2.2,2.3,1.9,2.3,1.7,
                1.5,1.6,1.5,2.1,1.6,1.8,2.1,1.7,1.5,1.6,1.5,2.0,2.1,1.7,2.1,1.5,
                2.2,2.3,2.2,2.8,2.3,2.5,2.8,2.4,2.2,2.3,2.2,2.7,2.8,2.4,2.8,2.2,
                2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
                1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
                2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
                1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
            ],
            // G C
            // C G
            [
                0.9,1.2,1.0,1.6,1.2,1.4,1.2,1.3,1.0,0.7,1.0,1.6,1.5,1.3,1.7,1.1,
                1.2,1.5,1.3,1.9,1.5,1.7,1.5,1.6,1.3,1.0,1.3,1.9,1.8,1.6,2.0,1.4,
                1.0,1.3,1.1,1.7,1.3,1.5,1.3,1.4,1.1,0.8,1.1,1.7,1.6,1.4,1.8,1.2,
                1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
                1.2,1.5,1.3,1.9,1.5,1.7,1.5,1.6,1.3,1.0,1.3,1.9,1.8,1.6,2.0,1.4,
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
                0.7,1.0,0.8,1.4,1.0,1.2,1.0,1.1,0.8,0.5,0.8,1.4,1.3,1.1,1.5,0.9,
                1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
                1.0,1.3,1.1,1.7,1.3,1.5,1.3,1.4,1.1,0.8,1.1,1.7,1.6,1.4,1.8,1.2,
                1.2,1.5,1.3,1.9,1.5,1.7,1.5,1.6,1.3,1.0,1.3,1.9,1.8,1.6,2.0,1.4,
                1.0,1.3,1.1,1.7,1.3,1.5,1.3,1.4,1.1,0.8,1.1,1.7,1.6,1.4,1.8,1.2,
                1.7,2.0,1.8,2.4,2.0,2.2,2.0,2.1,1.8,1.5,1.8,2.4,2.3,2.1,2.5,1.9,
                1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
                1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
                1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
                1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
            ],
            // G G
            // C C
            [
                1.0,1.3,1.1,1.8,1.3,1.5,0.9,1.3,1.1,1.2,1.1,1.6,1.6,1.3,1.5,1.1,
                1.3,1.6,1.4,2.1,1.6,1.8,1.2,1.6,1.4,1.5,1.4,1.9,1.9,1.6,1.8,1.4,
                1.1,1.4,1.2,1.9,1.4,1.6,1.0,1.4,1.2,1.3,1.2,1.7,1.7,1.4,1.6,1.2,
                1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
                1.3,1.6,1.4,2.1,1.6,1.8,1.2,1.6,1.4,1.5,1.4,1.9,1.9,1.6,1.8,1.4,
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
                0.8,1.1,0.9,1.6,1.1,1.3,0.7,1.1,0.9,1.0,0.9,1.4,1.4,1.1,1.3,0.9,
                1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
                1.1,1.4,1.2,1.9,1.4,1.6,1.0,1.4,1.2,1.3,1.2,1.7,1.7,1.4,1.6,1.2,
                1.3,1.6,1.4,2.1,1.6,1.8,1.2,1.6,1.4,1.5,1.4,1.9,1.9,1.6,1.8,1.4,
                1.1,1.4,1.2,1.9,1.4,1.6,1.0,1.4,1.2,1.3,1.2,1.7,1.7,1.4,1.6,1.2,
                1.8,2.1,1.9,2.6,2.1,2.3,1.7,2.1,1.9,2.0,1.9,2.4,2.4,2.1,2.3,1.9,
                1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
                1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
                1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
                1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
            ],
            // G T
            // C A
            [
                1.3,1.5,1.4,2.0,1.5,1.7,1.8,1.6,1.4,1.6,1.4,2.0,2.0,1.6,2.0,1.4,
                1.6,1.8,1.7,2.3,1.8,2.0,2.1,1.9,1.7,1.9,1.7,2.3,2.3,1.9,2.3,1.7,
                1.4,1.6,1.5,2.1,1.6,1.8,1.9,1.7,1.5,1.7,1.5,2.1,2.1,1.7,2.1,1.5,
                1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
                1.6,1.8,1.7,2.3,1.8,2.0,2.1,1.9,1.7,1.9,1.7,2.3,2.3,1.9,2.3,1.7,
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
                1.1,1.3,1.2,1.8,1.3,1.5,1.6,1.4,1.2,1.4,1.2,1.8,1.8,1.4,1.8,1.2,
                1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
                1.4,1.6,1.5,2.1,1.6,1.8,1.9,1.7,1.5,1.7,1.5,2.1,2.1,1.7,2.1,1.5,
                1.6,1.8,1.7,2.3,1.8,2.0,2.1,1.9,1.7,1.9,1.7,2.3,2.3,1.9,2.3,1.7,
                1.4,1.6,1.5,2.1,1.6,1.8,1.9,1.7,1.5,1.7,1.5,2.1,2.1,1.7,2.1,1.5,
                2.1,2.3,2.2,2.8,2.3,2.5,2.6,2.4,2.2,2.4,2.2,2.8,2.8,2.4,2.8,2.2,
                2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
                1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
                2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
                1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
            ],
            // T A
            // A T
            [
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
                2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
                2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
                2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
                2.2,2.3,2.2,2.8,2.3,2.5,2.8,2.4,2.2,2.3,2.2,2.7,2.8,2.4,2.8,2.2,
                2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
                2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
                2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
                2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
                2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
                2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
                2.4,2.5,2.4,3.0,2.5,2.7,3.0,2.6,2.4,2.5,2.4,2.9,3.0,2.6,3.0,2.4,
                1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
            ],
            // T C
            // A G
            [
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
                1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
                2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
                1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
                1.7,2.0,1.8,2.4,2.0,2.2,2.0,2.1,1.8,1.5,1.8,2.4,2.3,2.1,2.5,1.9,
                1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
                1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
                2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
                2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
                2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
                1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
                1.9,2.2,2.0,2.6,2.2,2.4,2.2,2.3,2.0,1.7,2.0,2.6,2.5,2.3,2.7,2.1,
                1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
            ],
            // T G
            // A C
            [
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
                1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
                2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
                1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
                1.8,2.1,1.9,2.6,2.1,2.3,1.7,2.1,1.9,2.0,1.9,2.4,2.4,2.1,2.3,1.9,
                1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
                1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
                2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
                2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
                2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
                1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
                2.0,2.3,2.1,2.8,2.3,2.5,1.9,2.3,2.1,2.2,2.1,2.6,2.6,2.3,2.5,2.1,
                1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
            ],
            // T T
            // A A
            [
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
                1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
                2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
                1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
                2.1,2.3,2.2,2.8,2.3,2.5,2.6,2.4,2.2,2.4,2.2,2.8,2.8,2.4,2.8,2.2,
                1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
                2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
                2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
                2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
                2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
                2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
                2.3,2.5,2.4,3.0,2.5,2.7,2.8,2.6,2.4,2.6,2.4,3.0,3.0,2.6,3.0,2.4,
                1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
            ],
    ],
        "22dh" :[
            // A A
            // T T
            [
                6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
                7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
                7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
                6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
                6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
                7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
                7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
                6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
                6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
                7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
                7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
                6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
                6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
                7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
                7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
                6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
            ],
            // A C
            // T G
            [
                -2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
                -1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
                -1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
                -1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
                -2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
                -1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
                -1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
                -1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
                -2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
                -1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
                -1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
                -1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
                -2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
                -1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
                -1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
                -1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
            ],
            // A G
            // T C
            [
                -0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
                0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
                0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
                -0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
                -0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
                0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
                0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
                -0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
                -0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
                0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
                0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
                -0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
                -0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
                0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
                0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
                -0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
            ],
            // A T
            // T A
            [
                8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
                9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
                8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
                8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
                8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
                9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
                8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
                8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
                8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
                9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
                8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
                8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
                8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
                9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
                8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
                8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
            ],
            // C A
            // G T
            [
                -2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
                -1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
                -2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
                -2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
                -2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
                -1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
                -2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
                -2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
                -2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
                -1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
                -2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
                -2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
                -2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
                -1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
                -2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
                -2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
            ],
            // C C
            // G G
            [
                -10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
                -10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
                -10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
                -11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
                -10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
                -10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
                -10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
                -11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
                -10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
                -10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
                -10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
                -11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
                -10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
                -10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
                -10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
                -11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
            ],
            // C G
            // G C
            [
                -9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
                -8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
                -9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
                -9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
                -9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
                -8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
                -9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
                -9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
                -9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
                -8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
                -9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
                -9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
                -9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
                -8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
                -9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
                -9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
            ],
            // C T
            // G A
            [
                -0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
                -0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
                -0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
                -1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
                -0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
                -0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
                -0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
                -1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
                -0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
                -0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
                -0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
                -1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
                -0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
                -0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
                -0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
                -1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
            ],
            // G A
            // C T
            [
                -3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
                -1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
                -2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
                1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
                -3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
                -1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
                -2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
                1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
                -3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
                -1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
                -2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
                1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
                -3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
                -1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
                -2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
                1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
            ],
            // G C
            // C G
            [
                -12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
                -9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
                -11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
                -7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
                -12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
                -9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
                -11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
                -7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
                -12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
                -9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
                -11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
                -7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
                -12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
                -9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
                -11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
                -7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
            ],
            // G G
            // C C
            [
                -10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
                -7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
                -9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
                -5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
                -10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
                -7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
                -9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
                -5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
                -10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
                -7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
                -9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
                -5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
                -10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
                -7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
                -9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
                -5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
            ],
            // G T
            // C A
            [
                -2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
                0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
                -1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
                3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
                -2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
                0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
                -1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
                3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
                -2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
                0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
                -1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
                3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
                -2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
                0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
                -1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
                3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
            ],
            // T A
            // A T
            [
                4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
                4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
                3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
                -3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
                4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
                4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
                3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
                -3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
                4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
                4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
                3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
                -3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
                4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
                4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
                3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
                -3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
            ],
            // T C
            // A G
            [
                -3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
                -3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
                -5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
                -12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
                -3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
                -3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
                -5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
                -12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
                -3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
                -3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
                -5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
                -12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
                -3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
                -3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
                -5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
                -12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
            ],
            // T G
            // A C
            [
                -2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
                -2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
                -3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
                -10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
                -2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
                -2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
                -3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
                -10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
                -2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
                -2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
                -3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
                -10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
                -2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
                -2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
                -3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
                -10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
            ],
            // T T
            // A A
            [
                6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
                6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
                4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
                -1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
                6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
                6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
                4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
                -1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
                6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
                6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
                4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
                -1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
                6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
                6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
                4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
                -1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
            ],
    ]
    
    
       


    }
    // Initialize variables
    let dH = 0;
    let dG = 0;

    let [seq1, seq2] = duplexStr.split("\n");
    let regionPosLi = loopRegionDict["region_pos"];
    let regionTypeLi = loopRegionDict["region_type"];
    let regionIdx = 0;

    while (regionIdx <= regionPosLi.length) {
        let start = regionIdx !== 0 ? regionPosLi[regionIdx - 1] : -1;
        let end = regionIdx !== regionPosLi.length ? regionPosLi[regionIdx] : seq1.length;

        // Stack regions
        if (regionIdx % 2 === 0 || regionIdx === regionPosLi.length) {
            start += 1;
            end -= 1;
            let stackLength = end - start + 1;
            if (stackLength > 1) {
                let segment = seq1.slice(start, end + 1);
                let [stackDh, stackDg] = stackEnergy(segment);
                dH += stackDh;
                dG += stackDg;
            }
        } else {
            let [iiDh, iiDg] = intermolecularInitiationEnergy();
            dH += iiDh;
            dG += iiDg;
            let loopIdx = Math.floor(regionIdx / 2);
            let loopType = regionTypeLi[loopIdx];
            let loopLength = end - start + 1;
            let segment1 = seq1.slice(start - 1, end + 2);
            let segment2 = seq2.slice(start - 1, end + 2);
            let segment = segment1.includes("-") ? segment2 : segment1;
            let seq = segment1.includes("-") ? reverseStr(seq2) : seq1;

            if (!loopType) {
                let [bulgeDh, bulgeDg] = bulgeEnergy(segment, loopLength, seq, start);
                dH += bulgeDh;
                dG += bulgeDg;
            } else {
                let [intLoopDh, intLoopDg] = intLoopEnergy(segment1, segment2, intLoopEnergyDict);
                dH += intLoopDh;
                dG += intLoopDg;
            }
        }
        regionIdx += 1;
    }

    // Calculate Tm
    console.log(`\ndG is ${dG}`);
    let dS = ((dH - dG) / 310.15) * 1000;
    dH *= 1000;
    let Tm = 0;
    console.log(`dH is ${dH}`);
    console.log(`dS is ${dS}`);

    // Initialize constants
    const T_KELVIN = 273.15;
    let K_mM = 50;
    let base = 4000000000;

    // Salt parameters
    const DNA_nM = 50;
    const dmsoConc = 0;
    const dmsoFact = 0.6;
    const formamideConc = 0.8;
    const divalent = 1.5;
    const dntp = 0.6;

    // Symmetry correction
    const sym = symmetry(seq1);
    if (sym) {
        dS += -1.4;
        base /= 4;
    }

    // Terminal AT penalty
    for (let i of [seq1[0], seq1[seq1.length - 1]]) {
        if (i === "A" || i === "T") {
            dS += 4.1;
            dH += 2300;
        } else {
            dS += -2.8;
            dH += 100;
        }
    }

    // GC content
    const GCCount = (countBase(seq1, "C") + countBase(seq1, "G") + countBase(seq2, "C") + countBase(seq2, "G")) / 2;
    K_mM += divalentToMonovalent(divalent, dntp);
    dS += 0.368 * (seq1.length - 1) * Math.log(K_mM / 1000);

    // Calculate Tm
    Tm = dH / (dS + 1.987 * Math.log(DNA_nM / base)) - T_KELVIN;
    Tm -= dmsoConc * dmsoFact;
    Tm += (0.453 * GCCount / seq1.length - 2.88) * formamideConc;
    console.log(`Tm is ${Tm}`);

    return Tm;
}

function main(){
    let duplex_str = "GCCCG\nCGG-C"
    let loop_region_dict = loopDetective(duplex_str)
    calcTmByNN(duplex_str, loop_region_dict)
}

main()