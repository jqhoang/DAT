package catanddog;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


/**
 * 
 * SequenceComparator is a local FASTA format file comparator, using REF as an indicator
 * for when the DNA sequence starts, two local files can be compared for percent identity
 * or similiarity. Does not account for gaps. This class is also capable of computing the RNA
 * amino acid composition.
 *
 * @author John A01001941
 * @version 2018
 */

public class SequenceComparator {
    
    private static File file1;
    private static File file2;
    public static Scanner scan;
    public static String[] aminoacids = {"Phenyl-alanine", "Leucine", "Serine", "Tyrosine", "Cysteine", 
            "Tryptophan", "Proline", "Histidine", "Glutamine", "Arginine", "Isoleucine",
            "Threonine", "Asparagine", "Lysine", "Valine", "Alanine", "Aspartic Acid", "Glutamic Acid",
            "Glycine"};
    public Scanner scanFile;
    public Scanner scanFile2;
    
    public SequenceComparator (File a, File b){
        
        file1 = a;
        file2 = b;
        
    }
    
    public static void displayMsg(String msg) {
        System.out.println(msg);
    }
    
    //does not account for N's
    public static String aminos(File file1) throws FileNotFoundException {
        
        /*int phe, leu, ser, tyr, cys, trp, pro, his, gin,
        arg, ile, met, thr, asn, lys, val, ala, asp, glu, gly = 0;*/
        String codon;
        String fullSeq = sequence(file1);
        if (fullSeq.length() % 3 != 0) {
            return "Sequence incompatible, indivisible by 3.";
        }
        String[] sequences = fullSeq.split("AUG");
        String[][] aainfo = new String[sequences.length][sequences.length];
        System.out.println(sequences.length);
        for (int k = 0; k < sequences.length; k++) {
            int[] aas = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            for (int i = 0; i < sequences[k].length();) {
                codon = sequences[k].substring(i, i + 1) + sequences[k].substring(i + 1, i + 2) + sequences[k].substring(i + 2, i + 3);
                if (codon.equals("UUU") || codon.equals("UUC")){
                    aas[0]++;
                    //phenyl-alanine
                } else if (codon.equals("UUA") || codon.equals("UUG") || (codon.equals("CUU") 
                        || codon.equals("CUC") || codon.equals("CUA") || codon.equals("CUG"))) {
                    aas[1]++;
                    //leucine
                } else if (codon.equals("UCU") || codon.equals("UCA") || codon.equals("UCC") 
                        || codon.equals("UCG") || codon.equals("AGU") || codon.equals("AGC")) {
                    aas[2]++;
                    //serine
                } else if (codon.equals("UAU") || codon.equals("UAC")) {
                    aas[3]++;
                    //tyrosine
                } else if (codon.equals("UGU") || codon.equals("UGC")) {
                    aas[4]++;
                    //cysteine
                } else if (codon.equals("UGG")) {
                    aas[5]++;
                    //tryptophan
                } else if (codon.equals("CCU") || codon.equals("CCC") || codon.equals("CCA") 
                        || codon.equals("CCG")) {
                    aas[6]++;
                    //proline
                } else if (codon.equals("CAU") || codon.equals("CAC")) {
                    aas[7]++;
                    //histidine
                } else if (codon.equals("CAA") || codon.equals("CAG")) {
                    aas[8]++;
                    //glutamine
                } else if (codon.equals("CGU") || codon.equals("CGC") || codon.equals("CGA") 
                        || codon.equals("CGG") || codon.equals("AGA") || codon.equals("AGG")) {
                    aas[9]++;
                    //arginine
                } else if (codon.equals("AUU") || codon.equals("AUC") || codon.equals("AUA")) {
                    aas[10]++;
                    //isoleucine
                } else if (codon.equals("ACU") || codon.equals("ACC") || codon.equals("ACA") 
                        || codon.equals("ACG")) {
                    aas[11]++;
                    //threonine
                } else if (codon.equals("AAU") || codon.equals("AAC")) {
                    aas[12]++;
                    //asparagine
                } else if (codon.equals("AAA") || codon.equals("AAG")) {
                    aas[13]++;
                    //lysine
                } else if (codon.equals("GUU") || codon.equals("GUC") || codon.equals("GUA") || codon.equals("GUG")) {
                    aas[14]++;
                    //valine
                } else if (codon.equals("GCU") || codon.equals("GCA") || codon.equals("GCC") || codon.equals("GCG")) {
                    aas[15]++;
                    //alanine
                } else if (codon.equals("GAU") || codon.equals("GAC")) {
                    aas[16]++;
                    //aspartic acid
                } else if (codon.equals("GAA") || codon.equals("GAG")) {
                    aas[17]++;
                    //glutamic acid
                } else if (codon.equals("GGU") || codon.equals("GGC") || codon.equals("GGA") || codon.equals("GGG")) {
                    aas[18]++;
                    //glycine
                } 
                i = i + 3;
            }
            String seqDesc = "";
            for (int l = 0; l < aas.length; l++) {
                seqDesc += " " + aminoacids[l] + " : " + aas[l];
            }
            aainfo[k][k] = sequences[k] + seqDesc;
            System.out.println(aainfo[k][k]);
        }
        return "holder";
    }
    
    private static String sequence(File file1) throws FileNotFoundException {
        
        String f1 = "";
        Scanner scanFile = new Scanner (file1);
        
        String REF1;
        
        REF1 = scanFile.nextLine();
        while (scanFile.hasNext()) {
            f1 += scanFile.next();
        }
        
        scanFile.close();
        return f1;

    }
    
    public static double percentIdentity(String file1, String file2) throws FileNotFoundException {

        double matches;
        double lengthOfRegion;

        String f1 = file1;
        String f2 = file2;
        int smallerString;
        int startPos;
        
        matches = 0;
        
        
        if (f1.length() > f2.length()) {
            smallerString = 2;
            lengthOfRegion = f2.length();
        } else {
            smallerString = 1;
            lengthOfRegion = f1.length();
        }
        
        startPos = 0;
        //one contains the other
        char[] convertedf1;
        char[] convertedf2;
        convertedf1 = new char[f1.length()];
        convertedf2 = new char[f2.length()];
        String convertedf1str = "";
        String convertedf2str = "";
        
        for (int i = 0; i < convertedf1.length; i++) {
            convertedf1[i] = f1.charAt(i);
        }
        for (int i = 0; i < convertedf2.length; i++) {
            convertedf2[i] = f2.charAt(i);
        }
        
        if (convertedf1.length <= convertedf2.length) {
            for (int i = 0; i < convertedf1.length; i++) {
                if (convertedf1[i] == 'N') {
                    convertedf1[i] = convertedf2[i];
                }
            }
        }
        if (convertedf2.length <= convertedf1.length) {
            for (int i = 0; i < convertedf2.length; i++) {
                if (convertedf2[i] == 'N') {
                    convertedf2[i] = convertedf1[i];
                }
            }
        }
        for (int k = 0; k < convertedf1.length; k++) {
            convertedf1str += convertedf1[k];
        }
        System.out.println(convertedf1str);

        for (int k = 0; k < convertedf2.length; k++) {
            convertedf2str += convertedf2[k];
        }
        System.out.println(convertedf2str);
        
        if (convertedf2str.contains(f1) || convertedf1str.contains(f2)) {
            return 100.00;
        }
        //making the string equivalents of the array
        
        // exact match
        if (smallerString == 1) {
            for (startPos = 0; startPos < convertedf1str.length(); startPos ++) {
                if (convertedf1str.charAt(startPos) == convertedf2str.charAt(startPos)) {
                    matches++;
                }
            }
        } else {
            for (startPos = 0; startPos < convertedf2str.length(); startPos ++) {
                if (convertedf2str.charAt(startPos) == convertedf1str.charAt(startPos)) {
                    matches++;
                }
            }
        }
        
        return matches/lengthOfRegion * 100;
        
    }
  


public static void main(String[] args) throws FileNotFoundException {
    scan = new Scanner(System.in);
    SequenceComparator.displayMsg("The first file path is:");
    file1 = new File(scan.nextLine());
    System.out.println(file1);
    SequenceComparator.displayMsg("The second file path is:");
    file2 = new File(scan.nextLine());
    System.out.println(file2);
    System.out.println(SequenceComparator.percentIdentity(SequenceComparator.sequence(file1), SequenceComparator.sequence(file2)));
    System.out.println("The file to grab the amino acid composition of is: ");
    File file3 = new File(scan.nextLine());
    aminos(file3);
    
    
}  
}