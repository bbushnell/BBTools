package prok;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import dna.AminoAcid;
import dna.GeneticCode;
import shared.Shared;
import stream.Read;
import structures.IntList;

/** Finite calling/translation/training oracles, independent of PGM biological accuracy.
 * @author Keqing
 */
public final class CallGenesGeneticCodeTest {

	public static void main(String[] args) throws Exception{
		Shared.setThreads(1);
		if(args.length!=1){throw new IllegalArgumentException("Expected fixture output directory");}
		final Path dir=Paths.get(args[0]);
		Files.createDirectories(dir);
		final byte[] legacy=AminoAcid.codeToByte.clone();
		final boolean[] starts=GeneModel.isStartCodon.clone(), stops=GeneModel.isStopCodon.clone();
		final GeneticCode eleven=GeneticCode.forTable(11), four=GeneticCode.forTable(4), twentyFive=GeneticCode.forTable(25);
		final Path customFile=dir.resolve("custom.tsv");
		Files.write(customFile, customTable().getBytes(StandardCharsets.US_ASCII));
		final GeneticCode custom=GeneticCode.load(customFile.toString());
		check(CallGenes.selectGeneticCode(null, null)==null, "No selector must preserve legacy policy");
		reject("mutually exclusive", ()->CallGenes.selectGeneticCode(11, customFile.toString()));
		reject("Unsupported translation table", ()->CallGenes.selectGeneticCode(999, null));
		check(CallGenes.selectGeneticCode(null, customFile.toString()).id()==0, "Custom selector lost custom identity");

		for(int strand=0; strand<2; strand++){
			checkOrf("ATGAAATGACCCTAA", null, strand, 0, 8, true, false, "MK");
			checkOrf("ATGAAATGACCCTAA", eleven, strand, 0, 8, true, false, "MK");
			checkOrf("ATGAAATGACCCTAA", four, strand, 0, 14, true, false, "MKWP");
			checkOrf("ATGAAATGACCCTAA", twentyFive, strand, 0, 14, true, false, "MKGP");
			checkOrf("TAAGTGCCCTAA", eleven, strand, 3, 11, false, false, "MP");
			checkOrf("TAAGTGCCCTAA", null, strand, 3, 11, false, false, "VP");
			checkOrf("GTGCCCTAA", eleven, strand, 0, 8, true, false, "VP");
			checkOrf("ATGAAACCC", four, strand, 0, 8, true, true, "MKP");
			checkOrf("ATGAAACCCT", four, strand, 0, 8, true, true, "MKP");
			checkOrf("ATGAAANNN", four, strand, 0, 8, true, true, "MKX");
			checkOrf("TAAAAACCCTAA", custom, strand, 3, 11, false, false, "MP");
			checkOrf("AAACCCTAA", custom, strand, 0, 8, true, false, "KP");
		}
		check(GeneCaller.makeOrfsForFrame("noStart", bytes("TAAAAACCCTAA"), 0, 0, 3).isEmpty(),
			"Legacy policy must not acquire the custom AAA start");
		check(!GeneCaller.validSyntheticEdgeStart(-64, four), "Ambiguous rolling codon must not enter strict code API");
		check(GeneCaller.validSyntheticEdgeStart(GeneticCode.codon("TGA"), four)
			&& !GeneCaller.validSyntheticEdgeStart(GeneticCode.codon("TGA"), eleven), "Edge eligibility must use selected stops");
		checkTraining(null, "findStartCodons", "CCCATTCCCATG", new int[]{9});
		checkTraining(eleven, "findStartCodons", "CCCATTCCCATG", new int[]{3,9});
		checkTraining(eleven, "findStopCodons", "CCCTGACCCTAA", new int[]{5,11});
		checkTraining(four, "findStopCodons", "CCCTGACCCTAA", new int[]{11});
		final byte[] wrongStop=bytes("ATGAAATGA");
		final Orf invalid=new Orf("invalid",0,8,1,0,wrongStop,true,ProkObject.CDS,false,false);
		AminoAcid.reverseComplementBasesInPlace(wrongStop);
		final byte[] before=wrongStop.clone();
		final Read invalidRead=new Read(wrongStop,null,"invalid",0);
		final ArrayList<Orf> invalidList=new ArrayList<Orf>();
		invalidList.add(invalid);
		reject("Complete CDS stop", ()->CallGenes.translate(invalidRead,invalidList,four));
		check(Arrays.equals(before,invalidRead.bases) && invalid.flipped()==1,
			"A rejected minus-strand translation must restore the input and ORF orientation");
		check(Arrays.equals(legacy, AminoAcid.codeToByte) && Arrays.equals(starts, GeneModel.isStartCodon)
			&& Arrays.equals(stops, GeneModel.isStopCodon), "Independent callers must not mutate global code tables");
		System.out.println("PASS CallGenesGeneticCodeTest: boundaries, strands, initiation, edges, custom code, training sites and isolation");
	}

	/** Coordinates and protein strings are hand-specified; source input must survive both-strand translation. */
	private static void checkOrf(String sequence, GeneticCode code, int strand, int start, int stop,
			boolean startTruncated, boolean stopTruncated, String protein){
		final byte[] biological=bytes(sequence);
		final ArrayList<Orf> found=GeneCaller.makeOrfsForFrame("fixture", biological, 0, strand, 3, code);
		check(found.size()==1, "Expected one frame0 ORF for "+sequence+", got "+found.size());
		final Orf orf=found.get(0);
		final int genomicStart=strand==0 ? start : biological.length-1-stop;
		final int genomicStop=strand==0 ? stop : biological.length-1-start;
		check(orf.start==genomicStart && orf.stop==genomicStop, "Wrong selected-code ORF coordinates: "+orf);
		check(orf.startTruncated==startTruncated && orf.stopTruncated==stopTruncated,
			"Biological boundary provenance changed for "+sequence);
		final byte[] genomic=biological.clone();
		if(strand==1){AminoAcid.reverseComplementBasesInPlace(genomic);}
		final Read input=new Read(genomic, null, "fixture", 0);
		final int flipped=orf.flipped();
		final ArrayList<Read> proteins=CallGenes.translate(input, found, code);
		check(proteins!=null && proteins.size()==1 && protein.equals(new String(proteins.get(0).bases, StandardCharsets.US_ASCII)),
			"Wrong protein for "+sequence+" under "+(code==null ? "legacy" : code.id()));
		check(orf.start==genomicStart && orf.stop==genomicStop && orf.flipped()==flipped,
			"Translation must restore ORF coordinates and flip state");
		final byte[] expected=biological.clone();
		if(strand==1){AminoAcid.reverseComplementBasesInPlace(expected);}
		check(Arrays.equals(expected, input.bases), "Translation must restore the genomic input sequence");
	}

	private static void checkTraining(GeneticCode code, String name, String sequence, int[] expected) throws Exception{
		final GeneModel model=new GeneModel(true, code);
		final Method method=GeneModel.class.getDeclaredMethod(name, byte[].class, IntList.class, BitSet.class);
		method.setAccessible(true);
		final IntList sites=new IntList();
		method.invoke(model, bytes(sequence), sites, new BitSet());
		check(sites.size==expected.length, "Wrong training site count for "+name);
		for(int i=0; i<expected.length; i++){check(sites.get(i)==expected[i], "Wrong training site for "+name+" at "+i);}
	}

	/** Complete synthetic code: AAA can initiate; internal AAA remains K; TGA elongates W. */
	private static String customTable(){
		final StringBuilder text=new StringBuilder("codon\tamino_acid\tstart\n");
		for(int i=0; i<64; i++){
			final String codon=AminoAcid.codonToString(i);
			text.append(codon).append('\t').append(codon.equals("TGA") ? 'W' : AminoAcid.toChar(codon)).append('\t')
				.append(codon.equals("ATG") || codon.equals("AAA") ? '1' : '0').append('\n');
		}
		return text.toString();
	}
	private static byte[] bytes(String value){return value.getBytes(StandardCharsets.US_ASCII);}
	private static void reject(String message, Runnable action){
		try{action.run();}catch(IllegalArgumentException e){check(e.getMessage().contains(message), "Wrong diagnostic: "+e); return;}
		throw new AssertionError("Expected diagnostic: "+message);
	}
	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}
}
