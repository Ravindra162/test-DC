
// Avalanche effect and Entropy Calculation......


const crypto = require('crypto');

// Function to calculate Hamming distance between two strings
function hammingDistance(str1, str2) {
    let distance = 0;
    const length = Math.min(str1.length, str2.length);
    for (let i = 0; i < length; i++) {
        if (str1[i] !== str2[i]) {
            distance++;
        }
    }
    distance += Math.abs(str1.length - str2.length); // Add the difference in lengths
    return distance;
}

// Function to calculate entropy of a string
function calculateEntropy(str) {
    const freqMap = {};
    for (let char of str) {
        freqMap[char] = (freqMap[char] || 0) + 1;
    }
    const length = str.length;
    let entropy = 0;
    for (let char in freqMap) {
        const probability = freqMap[char] / length;
        entropy -= probability * Math.log2(probability);
    }
    return entropy;
}

// Example usage with provided inputs
const originalMAC = "D0A1F5C9E7B4"; // Placeholder for original MAC
const RSencodedMAC = "AC0F6F82E5A9C3B7";  
let HashEncodedMAC;

// Measure avalanche effect (Hamming distance)
const hashOriginalMAC = crypto.createHash('sha256').update(originalMAC).digest('hex');
console.log("Hashed Original MAC:", hashOriginalMAC);

// Measure avalanche effect (Hamming distance) and calculate entropy
const avalancheEffect_RS = hammingDistance(originalMAC, RSencodedMAC) / Math.max(originalMAC.length, RSencodedMAC.length);
console.log("Avalanche Effect (Hamming Distance) - RS: ", avalancheEffect_RS);

// Generate hash of original MAC and store it in HashEncodedMAC
HashEncodedMAC = crypto.createHash('sha256').update(originalMAC).digest('hex');

const avalancheEffect_Hash = hammingDistance(originalMAC, HashEncodedMAC) / Math.max(originalMAC.length, HashEncodedMAC.length);
console.log("Avalanche Effect (Hamming Distance) - Hash: ", avalancheEffect_Hash);

const originalEntropy = calculateEntropy(originalMAC);
const RSencodedEntropy = calculateEntropy(RSencodedMAC);
const hashEncodedEntropy = calculateEntropy(HashEncodedMAC);
console.log("Original MAC Entropy: ", originalEntropy);
console.log("RS Encoded MAC Entropy: ", RSencodedEntropy);
console.log("Hash Encoded MAC Entropy: ", hashEncodedEntropy);
