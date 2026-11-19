package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import idaligner.QuantumAligner;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlIdentity;
import prok.NcrnaDdlSavedData;
import prok.PrepareNcrnaDdlCohorts;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Native direction/denominator and reused-position controls; inconsistent
 * captured gene sequences must fail before one can determine an identity bin.
 * @author Brian, Raiden
 */
public final class NcrnaDdlIdentityTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlIdentityTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown identity fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_fixture_directory without external inputs or sampling");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Identity fixture invariants require assertions");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Identity fixtures cannot overwrite existing evidence");
			alignment();
			conflictingGenes(root.getPath());
			coherentGenes(root.getPath(), false);
			coherentGenes(root.getPath(), true);
			pp.outstream.println("NCRNA_DDL_IDENTITY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void alignment(){
		final byte[] a=bytes("AAAAAAAAAAAAAAAAAAAA");
		final byte[] mismatch=bytes("AAAAAAAAACAAAAAAAAAA");
		final int[] scratch={987654, -9, 1234, -5};
		check(NcrnaDdlIdentity.align(a, a, scratch)==1f, "Identical complete genes have native consensus identity one");
		check(NcrnaDdlIdentity.align(a, mismatch, scratch)==19f/20,
			"One substitution in twenty complete consensus bases uses the native M/(M+S+I+D) denominator");
		final byte[] longer=bytes("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		final int[] direct=new int[4];
		direct[0]=Integer.MIN_VALUE;
		final float expected=QuantumAligner.alignStatic(longer, a, direct);
		Arrays.fill(scratch, Integer.MAX_VALUE);
		check(NcrnaDdlIdentity.align(longer, a, scratch)==expected && Arrays.equals(scratch, direct),
			"A longer consensus must remain the query; dirty endpoint scratch cannot turn into the next input offset");
		check(NcrnaDdlIdentity.align(a, a, scratch)==1f,
			"A short alignment after a longer one cannot reuse its bandwidth hint or returned endpoints");
	}

	private static void conflictingGenes(final String root){
		final RankSavedNcrnaWindowsTest.Fixture fixture=new RankSavedNcrnaWindowsTest.Fixture(root, "conflicting_genes");
		final NcrnaDdlSavedData data=NcrnaDdlSavedData.load(fixture.dir+"/panels.tsv", fixture.dir+"/selection.tsv", 2, 3);
		// This old ranking fixture intentionally uses A and G for the same parent
		// to exercise different event ranks. It is not coherent gene evidence.
		try{
			NcrnaDdlIdentity.measure(data, new byte[][]{bytes("AAAAAAAAAAAAAAAAAAAA"), bytes("GGGGGGGGGGGGGGGGGGGG")});
		}catch(IllegalArgumentException invalid){
			check(invalid.getMessage().contains("captured genes for one parent must agree"),
				"The fixture must reject conflicting genes, not fail an unrelated saved metadata control");
			return;
		}
		throw new AssertionError("Conflicting complete genes cannot be selected by read order");
	}

	private static void coherentGenes(final String root, final boolean undefined){
		final RankSavedNcrnaWindowsTest.Fixture fixture=new RankSavedNcrnaWindowsTest.Fixture(root,
			undefined ? "undefined_genes" : "coherent_genes");
		// Keep the old ranking witnesses but make the gene evidence coherent.
		// The old producer's exact FASTA remains beside the rewritten test input.
		final String path=fixture.dir+"/alpha.windows.fa.gz";
		final String original=fixture.dir+"/alpha.original.windows.fa.gz";
		check(new File(path).renameTo(new File(original)), "Preserve the exact synthetic input before changing its gene witness");
		try(Input in=new Input(original); Output out=new Output(path)){
			boolean parentZero=false;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				if(line.length>0 && line[0]=='>'){
					final int event=Parse.parseInt(line, 3, line.length);
					parentZero=event<2 || event>=5;
				}else if(parentZero){Arrays.fill(line, undefined ? (byte)'N' : (byte)'A');}
				out.writer.println(line);
			}
		}
		final NcrnaDdlSavedData data=NcrnaDdlSavedData.load(fixture.dir+"/panels.tsv", fixture.dir+"/selection.tsv", 2, 3);
		final byte[][] references={bytes("aaaaaaaaaaaaaaaaaaaa"), bytes("GGGGGGGGGGGGGGGGGGGG")};
		final NcrnaDdlIdentity.Locus[] loci=NcrnaDdlIdentity.measure(data, references);
		check(loci.length==6 && loci[0].model()==0 && loci[0].modelEvent()==1
			&& loci[0].containingWindows()==7 && loci[0].definedWindows()==(undefined ? 0 : 7),
			"Frozen model selection uses the best native rank then earliest event; every containing window is accounted for");
		check(loci[0].status().equals(undefined ? NcrnaDdlIdentity.UNDEFINED_GENE : NcrnaDdlIdentity.MEASURED)
			&& (undefined ? Float.isNaN(loci[0].identity()) : loci[0].identity()==1f && loci[0].geneEvent()==0),
			"The earliest complete gene is independent of model-event choice; undefined sequence stays explicit NA");
		check(loci[5].model()==0 && loci[5].modelEvent()==0 && loci[5].geneEvent()==0 && loci[5].identity()==1f,
			"The second panel's same parent name must retain its own gene and model event");
		for(int i=1; i<5; i++){
			check(loci[i].status().equals(NcrnaDdlIdentity.NOT_NATIVE_RANKABLE) && Float.isNaN(loci[i].identity()),
				"Ambiguous, missing-window, no-model and native-zero parents stay outside measured native-common identities");
		}
		if(!undefined){reuse(fixture, data, references, loci);}
		references[0]=bytes("AAAAAAAAANAAAAAAAAAA");
		final NcrnaDdlIdentity.Locus[] missing=NcrnaDdlIdentity.measure(data, references);
		check(missing[5].status().equals(NcrnaDdlIdentity.UNDEFINED_CONSENSUS) && Float.isNaN(missing[5].identity()),
			"An undefined fixed consensus cannot silently be normalized into a measured identity");
	}

	/** An accepted smaller population supplies five old values; one new gene aligns. */
	private static void reuse(final RankSavedNcrnaWindowsTest.Fixture fixture, final NcrnaDdlSavedData data,
			final byte[][] references, final NcrnaDdlIdentity.Locus[] original){
		final String dir=fixture.dir;
		for(final String name:new String[]{"panels", "selection"}){
			try(Input in=new Input(dir+"/"+name+".tsv"); Output out=new Output(dir+"/old."+name+".tsv")){
				out.writer.println(in.reader.nextLine());
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					if(row.length>5 && row[0]!='b'){out.writer.println(row);}
				}
			}
		}
		final NcrnaDdlSavedData old=NcrnaDdlSavedData.load(dir+"/old.panels.tsv", dir+"/old.selection.tsv", 2, 2);
		final NcrnaDdlIdentity.Locus[] prior=NcrnaDdlIdentity.measure(old, references);
		final NcrnaDdlIdentity.Locus[] extended=NcrnaDdlIdentity.extend(data, references, old, prior);
		for(int i=0; i<extended.length; i++){
			check(extended[i].reused()==(i<5)
				&& Float.floatToRawIntBits(extended[i].identity())==Float.floatToRawIntBits(original[i].identity())
				&& extended[i].status().equals(original[i].status()),
				"An extension reuses identical source values, including missing values, and measures only the newly admitted parent");
		}
		prior[0]=original[5];
		boolean rejected=false;
		try{NcrnaDdlIdentity.extend(data, references, old, prior);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Exchanging a prior identity's model/gene-event witnesses cannot silently reuse the value");
		PrepareNcrnaDdlCohorts.main(new String[]{"panels="+dir+"/old.panels.tsv", "selection="+dir+"/old.selection.tsv",
			"ref="+dir+"/ref.fa", "expectedmodels=2", "out="+dir+"/prior", "t=1"});
		PrepareNcrnaDdlCohorts.main(new String[]{"panels="+dir+"/panels.tsv", "selection="+dir+"/selection.tsv",
			"ref="+dir+"/ref.fa", "expectedmodels=2", "out="+dir+"/extended", "previous="+dir+"/prior",
			"previouspanels="+dir+"/old.panels.tsv", "previousselection="+dir+"/old.selection.tsv", "previousgroups=2", "t=1"});
		try(Input in=new Input(dir+"/extended.preparation.tsv")){
			in.reader.nextLine();
			check(new String(in.reader.nextLine(), StandardCharsets.UTF_8).equals("3\t6\t5\t1\t0"),
				"The real preparation CLI must restore five values and measure one new identity without sketch work");
		}
	}

	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
}
