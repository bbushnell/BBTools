package test.prok;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import prok.NcrnaAnnotatedCorpus;
import prok.NcrnaFamilySurvey;
import prok.NcrnaTextIO.Input;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Independent sequence, inherited-flag, duplicate-order and genomic-context checks.
 * @author Brian Bushnell, Raiden */
public final class NcrnaAnnotatedCorpusTest{

	public static void main(String[] args){
		//TODO: Inherited fixture discards PreParser's selected-log owner; an explicit
		//outstream file can remain open after an in-process call. Launcher exit reclaims
		//it; the current controls use borrowed stderr. Separate annotate/report only.
		args=new PreParser(args, NcrnaAnnotatedCorpusTest.class, false).args;
		check(args.length==1, "Pass a fresh fixture directory");
		final File dir=new File(args[0]);check(dir.mkdir(), "Fixture directory must be new");
		final String[] names={"tid_11_fixture_GCF_000001_1", "tid_22_fixture_GCF_000002_1", "tid_33_fixture_GCF_000003_1"};
		write(new File(dir, names[0]+".fna.gz"), ">c\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n");
		write(new File(dir, names[1]+".fna.gz"), ">d\nTTTTNNNNAAAACCCC\n");
		write(new File(dir, names[0]+".gff.gz"), GFF);
		write(new File(dir, names[1]+".gff.gz"), "d\tPGAP\triboswitch\t5\t12\t.\t+\t.\tID=b;Dbxref=RFAM:RF00174\n");
		write(new File(dir, names[2]+".gff.gz"), "z\tPGAP\tCDS\t1\t10\t.\t+\t0\tID=cds\n");
		final ByteBuilder roster=new ByteBuilder().append("index\tdomain\tphylum\tgenome\tgff\tcompressed_bytes\tmtime_ms\tfeature_rows\tcds_rows\tretained_layers\tcounted_loci\texact_duplicate_components\tlayers_without_noncoding_branch\n");
		for(int i=0; i<names.length; i++){
			roster.append(i).append("\tbacteria\tP1\t").append(names[i]).tab().append(new File(dir, names[i]+".gff.gz").getAbsolutePath()).append("\t0\t0\t0\t0\t0\t0\t0\t0\n");
		}
		write(new File(dir, "roster.tsv"), roster.toString());
		for(int threads:new int[]{1, 3}){
			final String prefix=new File(dir, "t"+threads).getPath();
			NcrnaAnnotatedCorpus.main(new String[]{"roster="+new File(dir, "roster.tsv"), "family=RF00174", "expectedgenomes=3",
				"expectedloci=6", "flank=5", "t="+threads, "outprefix="+prefix});
			verify(prefix);
		}
		for(String suffix:new String[]{".all.fa.gz", ".full.fa.gz", ".contexts.fa.gz", ".loci.tsv.gz", ".genomes.tsv", ".summary.tsv"}){
			compare(new File(dir, "t1"+suffix), new File(dir, "t3"+suffix));
		}
		final String a="c\tPGAP\triboswitch\t5\t12\t.\t+\t.\tID=a;Dbxref=RFAM:RF00174\n";
		final String b="c\tPGAP\tgene\t5\t12\t.\t+\t.\tID=parent;partial=true;pseudo=true\n"
			+"c\tPGAP\triboswitch\t5\t12\t.\t+\t.\tID=b;Parent=parent;Dbxref=RFAM:RF00174\n";
		for(int order=0; order<2; order++){
			final File f=new File(dir, "duplicates"+order+".gff");write(f, order==0 ? a+b : b+a);
			final ArrayList<NcrnaFamilySurvey.Locus> loci=NcrnaFamilySurvey.readLoci(f.getPath(), new long[6]);
			check(loci.size()==1 && loci.get(0).partial && loci.get(0).pseudo,
				"Exact-coordinate duplicate flags must be conservative in either source row order");
			check(loci.get(0).sourceIds.equals("a,b") && loci.get(0).sourceLines.contains(","),
				"Duplicate source references must survive alongside the single physical locus");
			check(loci.get(0).earliestSourceLine==(order==0 ? 1 : 2),
				"The representative source line is the numeric minimum supporting RNA row, independent of hash iteration");
		}
		final String branchA="c\tPGAP\tncRNA\t1\t12\t.\t+\t.\tID=A;Dbxref=RFAM:RF00174;partial=true;pseudo=true\n";
		final String branchB="c\tPGAP\tncRNA\t1\t26\t.\t+\t.\tID=B;Dbxref=RFAM:RF00174\n";
		final String exons="c\tPGAP\texon\t1\t4\t.\t+\t.\tParent=A,B\n"
			+"c\tPGAP\texon\t8\t12\t.\t+\t.\tParent=A\n"
			+"c\tPGAP\texon\t22\t26\t.\t+\t.\tParent=B\n";
		for(int order=0; order<2; order++){
			final File f=new File(dir, "siblings"+order+".gff");write(f, order==0 ? branchA+branchB+exons : exons+branchB+branchA);
			final ArrayList<NcrnaFamilySurvey.Locus> loci=NcrnaFamilySurvey.readLoci(f.getPath(), new long[6]);
			check(loci.size()==2, "Shared exons preserve both RNA branches");boolean aFound=false, bFound=false;
			for(NcrnaFamilySurvey.Locus l:loci){
				if(l.sourceIds.equals("A")){aFound=l.partial && l.pseudo && l.spans.equals("1-4,8-12");}
				if(l.sourceIds.equals("B")){bFound=!l.partial && !l.pseudo && l.spans.equals("1-4,22-26");}
			}
			check(aFound && bFound, "A shared exon must not leak its other parent's partial/pseudogene status into the selected branch, in either order");
		}
		verifyContigGc(dir);
		verifyNativeGcConsumer();
		System.out.println("ANNOTATED_CORPUS_TEST_PASS raw_conservation parent_flags duplicate_order reverse multipart context edge ambiguous zero_genome parallel contig_gc default_off");
	}

	/** Literal whole-contig fractions differ from the extracted locus/window fractions. */
	private static void verifyContigGc(final File dir){
		final ByteBuilder rosterText=new ByteBuilder("index\tdomain\tphylum\tgenome\tgff\tcompressed_bytes\tmtime_ms\tfeature_rows\tcds_rows\tretained_layers\tcounted_loci\texact_duplicate_components\tlayers_without_noncoding_branch\n");
		for(int i=0; i<3; i++){
			final String name="tid_"+(44+i)+"_fixture_GCF_00000"+(4+i)+"_1";
			write(new File(dir, name+".fna.gz"), ">mixed\nggccaaaannnn\n>fraction\ngttaann\n>undefined\nnnrryy\n>allgc\ncgcg\n");
			write(new File(dir, name+".gff.gz"),
			"mixed\tPGAP\triboswitch\t1\t4\t.\t+\t.\tID=plus;Dbxref=RFAM:RF00174\n"
			+"mixed\tPGAP\triboswitch\t5\t8\t.\t-\t.\tID=minus;Dbxref=RFAM:RF00174\n"
			+"fraction\tPGAP\triboswitch\t2\t3\t.\t+\t.\tID=fraction;Dbxref=RFAM:RF00174\n"
			+"undefined\tPGAP\triboswitch\t1\t3\t.\t-\t.\tID=undefined;Dbxref=RFAM:RF00174\n"
				+"allgc\tPGAP\triboswitch\t2\t3\t.\t+\t.\tID=allgc;Dbxref=RFAM:RF00174\n");
			rosterText.append(i).append("\tbacteria\tP1\t").append(name).tab().append(new File(dir, name+".gff.gz").getAbsolutePath()).append("\t0\t0\t0\t0\t0\t0\t0\t0\n");
		}
		final String roster=new File(dir, "gc.roster.tsv").getPath();
		write(new File(roster), rosterText.toString());
		for(int flank:new int[]{0, 2}){
			String prior=null;
			for(int threads:new int[]{1, 3}){
				final String prefix=new File(dir, "gc_f"+flank+"_t"+threads).getPath();
				runGc(roster, prefix, flank, threads, "gccontig=t");
				final int[] rows={0};
				RrnaResourceIO.stream(prefix+".contexts.fa.gz", r->{
					final String contig=RrnaResourceIO.field(r.id, "contig");
					final float expected=contig.equals("mixed") ? .5f : contig.equals("fraction") ? .2f : contig.equals("undefined") ? 0f : 1f;
					check(contig.equals("mixed") || contig.equals("fraction") || contig.equals("undefined") || contig.equals("allgc"), "Only declared literal contigs may appear");
					check(Float.floatToRawIntBits(Float.parseFloat(RrnaResourceIO.field(r.id, "contig_gc")))==Float.floatToRawIntBits(expected),
						"GC must use complete contig A/C/G/T counts before cropping or strand reversal: "+r.id);
					rows[0]++;
				});
				check(rows[0]==15, "GC annotation must retain all fifteen literal loci across three genomes");
				if(prior!=null){compareOutputs(prior, prefix);}
				prior=prefix;
			}
			final String defaultPrefix=new File(dir, "gc_default_f"+flank).getPath();
			final String falsePrefix=new File(dir, "gc_false_f"+flank).getPath();
			runGc(roster, defaultPrefix, flank, 1, null);
			runGc(roster, falsePrefix, flank, 1, "gccontig=f");
			compareOutputs(defaultPrefix, falsePrefix);
			RrnaResourceIO.stream(defaultPrefix+".contexts.fa.gz", r->check(!r.id.contains("contig_gc="), "Default output must omit GC metadata"));
			for(String suffix:SUFFIXES){
				if(!suffix.equals(".contexts.fa.gz")){compare(new File(defaultPrefix+suffix), new File(prior+suffix));}
			}
		}
	}
	/** Runs the real public entry point with an independently declared tiny corpus. */
	private static void runGc(String roster, String prefix, int flank, int threads, String gc){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("roster="+roster, "outprefix="+prefix,
			"family=RF00174", "expectedgenomes=3", "expectedloci=15", "flank="+flank, "t="+threads));
		if(gc!=null){args.add(gc);}
		NcrnaAnnotatedCorpus.main(args.toArray(new String[0]));
	}
	private static void compareOutputs(String a, String b){
		for(String suffix:SUFFIXES){compare(new File(a+suffix), new File(b+suffix));}
	}
	/** Literal fractions exercise the unchanged caller arithmetic and actual header consumer. */
	private static void verifyNativeGcConsumer(){
		try{
			final Method gc=prok.GeneCaller.class.getDeclaredMethod("calcContigGC", byte[].class);
			final Method parse=prok.NcrnaBoundaryVectorGen.class.getDeclaredMethod("parseHeaderFloat", String.class, String.class);
			gc.setAccessible(true);parse.setAccessible(true);
			final byte[] sparse=new byte[10000];Arrays.fill(sparse, (byte)'a');sparse[0]='g';
			final byte[][] bases={new byte[0], "nry".getBytes(StandardCharsets.US_ASCII),
				"gcUuN".getBytes(StandardCharsets.US_ASCII), sparse};
			final float[] expected={0f, 0f, .5f, 1e-4f};
			for(int i=0; i<bases.length; i++){
				final byte[] original=bases[i].clone();
				final float actual=(Float)gc.invoke(null, bases[i]);
				check(Float.floatToRawIntBits(actual)==Float.floatToRawIntBits(expected[i]) && Arrays.equals(original, bases[i]),
					"Caller GC must preserve complete input bases and the literal defined-base convention; case="+i);
				final float restored=(Float)parse.invoke(null, "q contig_gc="+Float.toString(actual)+" trailing=kept", "contig_gc=");
				check(Float.floatToRawIntBits(restored)==Float.floatToRawIntBits(expected[i]), "Actual boundary consumer must preserve GC float bits including exponent form");
			}
			check(Float.isNaN((Float)parse.invoke(null, "q no_gc=true", "contig_gc=")), "Missing optional GC stays unavailable rather than becoming window GC");
			System.out.println("ANNOTATED_CORPUS_NATIVE_GC_PASS fractions=4 missing_field=1");
		}catch(ReflectiveOperationException failure){throw new IllegalStateException("Native GC consumer fixture failed", failure);}
	}

	private static void verify(final String prefix){
		final ObjectMap<String, String> exact=new ObjectMap<String, String>(String.class, String.class);
		RrnaResourceIO.stream(prefix+".all.fa.gz", r->{check(exact.put(RrnaResourceIO.first(r.id), new String(r.bases, StandardCharsets.US_ASCII))==null, "Extracted locus IDs must be unique");});
		final LineParser1 p=new LineParser1('\t');
		int rows=0, full=0;boolean minus=false, multipart=false, edge=false, partial=false, pseudo=false, ambiguous=false;
		try(Input owner=new Input(prefix+".loci.tsv.gz")){
			final ByteFile in=owner.reader;
			check(in.nextLine()!=null, "Corpus ledger needs a schema header");
			for(byte[] row; (row=in.nextLine())!=null;){
				p.set(row);check(p.terms()==21, "Corpus ledger must preserve provenance and flag fields");rows++;
				final String seq=exact.get(p.parseString(0)), spans=p.parseString(8), status=p.parseString(14);
				if(status.equals("ANNOTATION_FULL_COPY")){full++;}
				if(spans.equals("33-40")){minus=seq.equals("GGGGTTTT") && p.termEquals("-", 7) && status.equals("ANNOTATION_FULL_COPY");}
				if(spans.equals("49-52,57-60")){multipart=seq.equals("CCCCTTTT") && status.equals("MULTIPART_UNRESOLVED");}
				if(spans.equals("1-3")){edge=seq.equals("AAA") && p.termEquals("true", 16) && p.parseInt(18)==0 && p.parseInt(19)==5;}
				if(spans.equals("21-28")){partial=p.termEquals("true", 12) && status.equals("SOURCE_PARTIAL");}
				if(spans.equals("5-12") && p.termEquals("c", 6)){pseudo=p.termEquals("true", 13) && status.equals("PSEUDOGENE") && seq.equals("CCCCGGGG");}
				if(p.termEquals("d", 6)){ambiguous=seq.equals("NNNNAAAA") && p.parseInt(17)==4 && status.equals("ANNOTATION_FULL_COPY");}
			}
		}
		check(rows==6 && exact.size()==6 && full==3 && minus && multipart && edge && partial && pseudo && ambiguous,
			"Independent literal sequences and source-flag dispositions must all match; siblings must not inherit each other's flags");
		final int[] fullCount={0};RrnaResourceIO.stream(prefix+".full.fa.gz", r->{fullCount[0]++;});
		check(fullCount[0]==3, "Only the three annotation-full simple loci enter full.fa");
		final int[] contextCount={0};RrnaResourceIO.stream(prefix+".contexts.fa.gz", r->{contextCount[0]++;
			if(r.id.contains("start=49 stop=60")){
				check(r.bases.length==21 && r.id.contains("lflank=4 rflank=5"), "Multipart context remains contiguous, including gap and clipped reverse flanks");
			}
		});check(contextCount[0]==6, "Contexts retain every raw physical locus");
	}
	private static void write(final File path, final String text){
		//TODO: Inherited helper starts before ownership and closes only after a successful print.
		//Annotate/report under the current scope; new GC cases did not change this output helper.
		final ByteStreamWriter out=RrnaResourceIO.writer(path.getPath());out.print(text);RrnaResourceIO.finish(out);
	}
	private static void compare(final File a, final File b){
		check(Arrays.equals(decoded(a), decoded(b)), "Exact decompressed outputs must match: "+a+" / "+b);
	}
	/** Native raw input keeps terminators in the parity comparison; fixtures are tiny. */
	private static byte[] decoded(File file){
		final ByteBuilder bytes=new ByteBuilder();final byte[] buffer=new byte[4096];
		try(InputStream in=ReadWrite.getInputStream(file.getPath(), true, false)){
			for(int n; (n=in.read(buffer))>=0;){if(n>0){bytes.append(buffer, 0, n);}}
		}catch(IOException failure){throw new IllegalStateException("Fixture read failed: "+file, failure);}
		return bytes.toBytes();
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
	private static final String[] SUFFIXES={".all.fa.gz", ".full.fa.gz", ".contexts.fa.gz", ".loci.tsv.gz", ".genomes.tsv", ".summary.tsv"};
	private static final String GFF=
		"c\tPGAP\triboswitch\t5\t12\t.\t+\t.\tID=a;Dbxref=RFAM:RF00174\n"
		+"c\tPGAP\triboswitch\t5\t12\t.\t+\t.\tID=alias;Dbxref=RFAM:RF00174;pseudo=true\n"
		+"c\tPGAP\tgene\t21\t28\t.\t+\t.\tID=p;partial=true\n"
		+"c\tPGAP\triboswitch\t21\t28\t.\t+\t.\tID=r;Parent=p;Dbxref=RFAM:RF00174\n"
		+"c\tPGAP\triboswitch\t33\t40\t.\t-\t.\tID=m;Dbxref=RFAM:RF00174\n"
		+"c\tPGAP\tncRNA\t49\t60\t.\t-\t.\tID=s;Dbxref=RFAM:RF00174\n"
		+"c\tPGAP\texon\t49\t52\t.\t-\t.\tParent=s\n"
		+"c\tPGAP\texon\t57\t60\t.\t-\t.\tParent=s\n"
		+"c\tPGAP\triboswitch\t1\t3\t.\t+\t.\tID=edge;Dbxref=RFAM:RF00174\n";
}
