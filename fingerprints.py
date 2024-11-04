"Class to handle Multiple Fingerprint"

from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
import pubchempy as pcp  # noqa
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
import numpy as np
import os

## JVM to run CDK
from jpype import isJVMStarted, startJVM, getDefaultJVMPath, JPackage

## Connecting with JAR file
if not isJVMStarted():
    cdk_path = os.path.join(os.path.dirname(__file__), "CDK", "cdk-2.9.jar")
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % cdk_path)
    cdk = JPackage("org").openscience.cdk
else:
    cdk = JPackage("org").openscience.cdk


class Fingerprint:
    def __init__(
        self,
    ):
        pass

    def calculate_fingerprint(self, smiles):
        pass

    @staticmethod
    def smiles_to_molecule(smiles: str):
        return AllChem.MolFromSmiles(smiles)


class RDKIT_MACCS(Fingerprint):
    def __init__(
        self,
    ):
        pass

    def calculate_fingerprint(self, smiles):
        """Calculate MACCS Fingerprint"""
        mol = RDKIT_MACCS.smiles_to_molecule(smiles)
        fp = MACCSkeys.GenMACCSKeys(mol)
        return fp


class PubChem(Fingerprint):
    def __init__(self):
        pass

    @staticmethod
    def calculate_fingerprint(smiles):
        """Calculate PubChem Fingerprint

        Source: https://github.com/deepchem/deepchem/blob/master/deepchem/feat/molecule_featurizers/pubchem_fingerprint.py
        """
        pubchem_compound = pcp.get_compounds(smiles, "smiles")
        if pubchem_compound:
            feature = [int(bit) for bit in pubchem_compound[0].cactvs_fingerprint]

            # Create an ExplicitBitVect of appropriate length
            fp = ExplicitBitVect(len(feature))
            for i, bit in enumerate(feature):
                if bit:
                    fp.SetBit(int(i))  # Explicitly cast i to int
            return fp
        else:
            return None  # Handle case when no compound is found


class ECFP6(Fingerprint):
    def __init__(
        self,
    ):
        pass

    @staticmethod
    def calculate_fingerprint(smiles):
        """Calculate ECFP6 Fingerprint"""
        molecule = ECFP6.smiles_to_molecule(smiles)
        return AllChem.GetMorganFingerprintAsBitVect(molecule, 3, nBits=1024)


class KlekotaRothFingerprint(Fingerprint):
    def __init__(self, nbits=4860):
        super().__init__()

        self.cdk_smile_parser = cdk.smiles.SmilesParser(
            cdk.DefaultChemObjectBuilder.getInstance()
        )
        self.cdk_klekota_roth_fingerprinter = cdk.fingerprint.KlekotaRothFingerprinter()
        self.nbits = nbits

    def smiles_to_molecule(self, smiles: str):
        mol = self.cdk_smile_parser.parseSmiles(smiles)
        return mol

    def calculate_fingerprint(self, smiles):
        """Calculate KlekotaRoth Fingerprint"""
        mol = self.smiles_to_molecule(smiles)

        # Create an instance of the fingerprinter
        fingerprinter = self.cdk_klekota_roth_fingerprinter

        fp = fingerprinter.getBitFingerprint(mol)

        # Create feature array directly from the set bits
        feature = np.zeros(self.nbits, dtype=int)
        bits = list(fp.getSetbits())

        feature[bits] = 1

        # Create and return an ExplicitBitVect using RDKit
        final_fp = ExplicitBitVect(self.nbits)
        for i in range(self.nbits):
            if feature[i]:
                final_fp.SetBit(i)

        return final_fp


class CDK_PUBCHEM(Fingerprint):
    def __init__(self, nbits=2048):
        self.cdk_smile_parser = cdk.smiles.SmilesParser(
            cdk.DefaultChemObjectBuilder.getInstance()
        )
        self.pubchem_fingerprinter = cdk.fingerprint.PubchemFingerprinter(
            cdk.silent.SilentChemObjectBuilder.getInstance()
        )
        self.nbits = nbits

    def smiles_to_molecule(self, smiles: str):
        mol = self.cdk_smile_parser.parseSmiles(smiles)
        return mol

    def calculate_fingerprint(self, smiles):
        """Calculate PubChem Fingerprint"""
        mol = self.smiles_to_molecule(smiles)

        # Create an instance of the fingerprinter
        fingerprinter = self.pubchem_fingerprinter

        fp = fingerprinter.getBitFingerprint(mol)

        # Create feature array directly from the set bits
        feature = np.zeros(self.nbits, dtype=int)
        bits = list(fp.getSetbits())

        feature[bits] = 1

        # Create and return an ExplicitBitVect using RDKit
        final_fp = ExplicitBitVect(self.nbits)
        for i in range(self.nbits):
            if feature[i]:
                final_fp.SetBit(i)

        return final_fp
