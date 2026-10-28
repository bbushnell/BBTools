package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.FamilyResourceResolver;
import prok.NcrnaFamily;
import prok.NcrnaModelThresholds;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.RrnaResourceIO;

/** Padded windows are not gene lengths; rejecting one model must preserve a
 * different passing model. Optional identity columns retain the same aligner.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaModelLengthBoundsTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");
		final Path root=Paths.get(args[0]);Files.createDirectory(root);final Path table=root.resolve("lengths.tsv");
		final byte[] gene=new byte[120], alphabet={'A','C','G','T'};final Random random=new Random(74063);
		for(int i=0; i<gene.length; i++){gene[i]=alphabet[random.nextInt(4)];}
		final byte[] high=gene.clone(), low=gene.clone();
		for(int i=40; i<120; i+=10){high[i]=NcrnaModelThresholdsTest.other(high[i]);}
		for(int i=40; i<120; i+=5){low[i]=NcrnaModelThresholdsTest.other(low[i]);}
		final byte[] bases=new byte[200];Arrays.fill(bases, (byte)'N');System.arraycopy(gene, 0, bases, 40, 120);
		final NcrnaFamily family=NcrnaModelThresholdsTest.family(high, low, gene);
		family.reuseConsensusAlignment=false;
		final ArrayList<Orf> baseline=NcrnaModelThresholdsTest.call(family, bases);
		check(baseline.size()==1 && baseline.get(0).trnaModel.equals("high"), "Unbounded generic verifier must select the higher-identity model");
		write(table, "high\t120\t140\t200\nlow\t120\t100\t130\n");
		NcrnaModelThresholds.loadLengths(table.toString(), family);
		check(!family.reuseConsensusAlignment && family.modelThresholds.sameIdentities(.6f, .6f), "Length metadata must not switch alignment or identity policy");
		final NcrnaScavenger bounded=RrnaEndpointCallerFeaturesTest.scavenger(family);bounded.setMappingOracleExhaustive(true);
		check(bounded.modelLengthFeasible(1, 1000), "A long padded window can contain a short valid gene");
		check(!bounded.modelLengthFeasible(0, 100), "A window shorter than the minimum admitted span can skip alignment");
		check(bounded.modelAlignedSpanAllowed(1, 100) && bounded.modelAlignedSpanAllowed(1, 130)
			&& !bounded.modelAlignedSpanAllowed(1, 99) && !bounded.modelAlignedSpanAllowed(1, 131), "Inclusive observed-span limits must reject only values outside the resource interval");
		final ArrayList<Orf> calls=bounded.scavenge("contig", bases, 0, new ArrayList<int[]>());
		check(calls.size()==1 && calls.get(0).trnaModel.equals("low") && bounded.modelSpanRejects()>0,
			"Higher-identity model with an invalid observed span must not hide a passing runner-up");
		final NcrnaModelThresholds good=family.modelThresholds;
		for(String invalid:new String[]{"high\t119\t100\t130\nlow\t120\t100\t130\n", "high\t120\t100\t130\n",
			"high\t120\t100\t130\nhigh\t120\t100\t130\n", "high\t120\t130\t100\nlow\t120\t100\t130\n"}){
			write(table, invalid);NcrnaModelThresholdsTest.reject(()->NcrnaModelThresholds.loadLengths(table.toString(), family));
			check(family.modelThresholds==good, "Invalid name/order/length/count data must not partially alter a live family");
		}
		final FamilyResourceResolver resolver=new FamilyResourceResolver();final FamilyResourceResolver.Overrides overrides=new FamilyResourceResolver.Overrides();
		overrides.root=root.resolve("absent_manifest").toString();
		check(resolver.optional("euk5S", "length_bounds", overrides)==null, "Optional lengths must not require a manifest for fully explicit legacy resource overrides");
		NcrnaModelThresholdsTest.reject(()->resolver.resolve("euk5S", "consensus", overrides));
		checkIdentities(table, high, low, gene, bases);
		System.out.println("NCRNA_MODEL_LENGTH_BOUNDS_TEST_PASS");
	}
	private static void checkIdentities(Path table, byte[] high, byte[] low, byte[] gene, byte[] bases){
		final NcrnaFamily family=NcrnaModelThresholdsTest.family(high, low, gene);
		family.reuseConsensusAlignment=false;
		final String plain="high\t120\t100\t130\nlow\t120\t100\t130\n";
		final String specific="high\t120\t100\t130\t1\nlow\t120\t100\t130\t.6\n";
		write(table, plain);NcrnaModelThresholds.loadLengths(table.toString(), family);
		final String baseline=RrnaEndpointCallerFeaturesTest.calls(NcrnaModelThresholdsTest.call(family, bases));
		writeIdentity(table, "high\t120\t100\t130\t.6\nlow\t120\t100\t130\t.6\n");
		NcrnaModelThresholds.loadLengths(table.toString(), family);
		check(baseline.equals(RrnaEndpointCallerFeaturesTest.calls(NcrnaModelThresholdsTest.call(family, bases))),
			"Explicit common cutoffs must preserve the same generic calls as absent identity columns");
		writeIdentity(table, specific);NcrnaModelThresholds.loadLengths(table.toString(), family);
		for(boolean quantum:new boolean[]{true, false}){
			for(boolean ranked:new boolean[]{true, false}){
				final NcrnaScavenger caller=RrnaEndpointCallerFeaturesTest.scavenger(family);
				caller.quantumThresh=quantum ? 120 : Integer.MAX_VALUE;caller.rankedModelFallback=ranked;
				final ArrayList<Orf> calls=caller.scavenge("contig", bases, 0, new ArrayList<int[]>());
				check(calls.size()==1 && calls.get(0).trnaModel.equals("low"),
					"Generic Quantum/Scrabble and first-passing policies must retain the eligible lower-identity model");
			}
		}
		final NcrnaModelThresholds good=family.modelThresholds;
		for(String bad:new String[]{"NaN", "-0.1", "1.1", ""}){
			writeIdentity(table, "high\t120\t100\t130\t"+bad+"\nlow\t120\t100\t130\t.6\n");
			if(bad.isEmpty()){
				boolean rejected=false;
				try{NcrnaModelThresholds.loadLengths(table.toString(), family);}
				catch(AssertionError expected){rejected=expected.getMessage().contains("numeric fields must be nonempty");}
				check(rejected, "savedDouble's explicit nonempty-field assertion must reject a blank identity");
			}else{NcrnaModelThresholdsTest.reject(()->NcrnaModelThresholds.loadLengths(table.toString(), family));}
			check(family.modelThresholds==good, "Malformed identity data must fail before replacing the live immutable profile");
		}
		writeIdentity(table, "high\t120\t100\t130\t.55\nlow\t120\t100\t130\t.6\n");
		NcrnaModelThresholds.loadLengths(table.toString(), family);
		check(family.modelThresholds.pass(0)==.55f && family.modelThresholds.borderline(0)==.55f,
			"A lower direct cutoff must not be hidden behind the inherited .6 rescue floor");
	}
	private static void writeIdentity(Path path, String rows){
		final ByteStreamWriter out=writer(path);
		try{out.print("model\tmodel_length\tmin_length\tmax_length\tidentity\n"+rows);}finally{RrnaResourceIO.finish(out);}
	}
	private static void write(Path path, String rows){
		final ByteStreamWriter out=writer(path);
		try{out.print("model\tmodel_length\tmin_length\tmax_length\n"+rows);}finally{RrnaResourceIO.finish(out);}
	}
	/** These files are deliberately replaced between cases inside the fresh test directory. */
	private static ByteStreamWriter writer(Path path){
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path.toString(), FileFormat.TEXT, null, true, true, false, false));
		out.start();return out;
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
