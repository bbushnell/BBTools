package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaAnnotatedCorpus;
import prok.NcrnaFamilySurvey;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Independent sequence, inherited-flag, duplicate-order and genomic-context checks.
 * @author Brian Bushnell, Raiden */
public final class NcrnaAnnotatedCorpusTest {

	public static void main(final String[] args){
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
		System.out.println("ANNOTATED_CORPUS_TEST_PASS raw_conservation parent_flags duplicate_order reverse multipart context edge ambiguous zero_genome parallel");
	}

	private static void verify(final String prefix){
		final ObjectMap<String, String> exact=new ObjectMap<String, String>(String.class, String.class);
		RrnaResourceIO.stream(prefix+".all.fa.gz", r->{check(exact.put(RrnaResourceIO.first(r.id), new String(r.bases, StandardCharsets.US_ASCII))==null, "Extracted locus IDs must be unique");});
		final ByteFile in=RrnaResourceIO.open(prefix+".loci.tsv.gz");final LineParser1 p=new LineParser1('\t');
		int rows=0, full=0;boolean minus=false, multipart=false, edge=false, partial=false, pseudo=false, ambiguous=false;
		try{
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
		}finally{check(!in.close(), "Corpus ledger read failed");}
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
		final ByteStreamWriter out=RrnaResourceIO.writer(path.getPath());out.print(text);RrnaResourceIO.finish(out);
	}
	private static void compare(final File a, final File b){
		final ByteFile x=RrnaResourceIO.open(a.getPath()), y=RrnaResourceIO.open(b.getPath());
		try{for(byte[] row; (row=x.nextLine())!=null;){check(Arrays.equals(row, y.nextLine()), "One/three-worker decoded output must match");}check(y.nextLine()==null, "Parallel output has additional rows");}
		finally{check(!x.close() && !y.close(), "Comparison input failed");}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
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
