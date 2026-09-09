package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;

/** Focused Java-8 fixture for the typed expected-copy producer. */
public final class MagQCExpectedCopyTableTypedTest {
	private MagQCExpectedCopyTableTypedTest(){}

	public static void main(String[] args) throws Exception{
		Path dir=Files.createTempDirectory("expected-copy-typed.");
		try{
			Paths p=inputs(dir);
			Path out=dir.resolve("table.tsv");
			MagQCExpectedCopyTableTyped.Table t=MagQCExpectedCopyTableTyped.build(p.counts.toString(),p.roster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString(),out.toString());
			check(t.population("global","-").denominator==3,"global includes partial and unknown");
			check(t.population("phylum","Firmicutes").denominator==1,"classified phylum only");
			check(close(t.expectation("global","-","P","0"),2.0),"global P0");
			check(close(t.expectation("global","-","P","1"),1.0/3.0),"global P1");
			check(close(t.expectation("global","-","N","r16"),1.0),"global r16");
			check(close(t.expectation("phylum","Firmicutes","P","1"),1.0),"phylum P1");
			check(close(t.expectation("global","-","P","0"),t.expectation("phylum","Firmicutes","P","0")),"typed lookup");
			MagQCExpectedCopyTableTyped.Table loaded=MagQCExpectedCopyTableTyped.load(out.toString(),p.counts.toString(),p.roster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString());
			check(close(loaded.expectation("global","-","N","trna"),1.0/3.0),"round trip");
			testInputOrder(dir,p,t);
			testHashMismatch(dir,p,out);
			testMalformed(dir,p,out);
			System.out.println("PASS MagQCExpectedCopyTableTypedTest");
		}finally{deleteTree(dir);}
	}

	private static void testInputOrder(Path dir,Paths p,MagQCExpectedCopyTableTyped.Table base) throws Exception{
		String ua=unit("tid_101_a.fa"),ub=unit("tid_102_b.fa"),uc=unit("tid_103_c.fa");
		Path shuffled=write(dir.resolve("counts-shuffled.tsv"),"#schema_version\tclean_typed_counts_v1\n#columns\t"+MagQCExpectedCopyTableTyped.COUNT_COLUMNS+"\n"
				+uc+"\tgenome_"+uc+"\t103\tN\ttrna\t0\n"
				+ua+"\tgenome_"+ua+"\t101\tP\t0\t2\n"
				+ua+"\tgenome_"+ua+"\t101\tP\t1\t1\n"
				+ub+"\tgenome_"+ub+"\t102\tN\tr16\t2\n"
				+ua+"\tgenome_"+ua+"\t101\tN\tr16\t1\n"
				+uc+"\tgenome_"+uc+"\t103\tN\tr16\t0\n"
				+ub+"\tgenome_"+ub+"\t102\tN\ttrna\t0\n"
				+ua+"\tgenome_"+ua+"\t101\tN\ttrna\t1\n"
				+ub+"\tgenome_"+ub+"\t102\tN\tr23\t0\n"
				+uc+"\tgenome_"+uc+"\t103\tN\tr23\t0\n"
				+ua+"\tgenome_"+ua+"\t101\tN\tr23\t1\n"
				+ub+"\tgenome_"+ub+"\t102\tN\tr5\t0\n"
				+uc+"\tgenome_"+uc+"\t103\tN\tr5\t0\n"
				+ua+"\tgenome_"+ua+"\t101\tN\tr5\t0\n"
				+ub+"\tgenome_"+ub+"\t102\tN\trother\t0\n"
				+uc+"\tgenome_"+uc+"\t103\tN\trother\t0\n"
				+ua+"\tgenome_"+ua+"\t101\tN\trother\t0\n"
				+ub+"\tgenome_"+ub+"\t102\tP\t0\t4\n");
		Path out=dir.resolve("shuffled.tsv");MagQCExpectedCopyTableTyped.Table t=MagQCExpectedCopyTableTyped.build(shuffled.toString(),p.roster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString(),out.toString());
		for(MagQCExpectedCopyTableTyped.Item item:base.items){check(close(t.expectation("global","-",item.type,item.key),base.expectation("global","-",item.type,item.key)),"input order global "+item.id());check(close(t.expectation("phylum","Firmicutes",item.type,item.key),base.expectation("phylum","Firmicutes",item.type,item.key)),"input order phylum "+item.id());}
	}

	private static void testHashMismatch(Path dir,Paths p,Path out) throws Exception{
		final String original=read(p.counts);
		write(p.counts,original+"\n");
		try{expectFailureContaining(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.load(out.toString(),p.counts.toString(),p.roster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString());}},"same-path count hash mismatch","provenance hash mismatch");}
		finally{write(p.counts,original);}
		final Path alternate=write(dir.resolve("counts-alternate.tsv"),original);
		expectFailureContaining(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.load(out.toString(),alternate.toString(),p.roster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString());}},"count path mismatch","provenance path mismatch");
	}

	private static void testMalformed(Path dir,Paths p,Path out) throws Exception{
		String counts=read(p.counts);
		String ua=unit("tid_101_a.fa"),ub=unit("tid_102_b.fa"),uc=unit("tid_103_c.fa");
		Path duplicate=write(dir.resolve("counts-duplicate.tsv"),counts+ua+"\tgenome_"+ua+"\t101\tN\tr16\t1\n");
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.build(duplicate.toString(),p.roster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString(),dir.resolve("bad-duplicate.tsv").toString());}},"duplicate typed count");
		Path missing=write(dir.resolve("counts-missing.tsv"),counts.replace(uc+"\tgenome_"+uc+"\t103\tN\ttrna\t0\n",""));
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.build(missing.toString(),p.roster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString(),dir.resolve("bad-missing.tsv").toString());}},"missing ncRNA count");
		Path badRoster=write(dir.resolve("roster-duplicate-tid.tsv"),read(p.roster).replace(ub+"\tgenome_"+ub+"\ttid_102_b.fa\t/assemblies/b.fa\t"+repeat('b')+"\t102", ub+"\tgenome_"+ub+"\ttid_102_b.fa\t/assemblies/b.fa\t"+repeat('b')+"\t101"));
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.build(p.counts.toString(),badRoster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString(),dir.resolve("bad-roster.tsv").toString());}},"duplicate roster tid");
		Path duplicateExcludedLabel=write(dir.resolve("labels-duplicate-excluded-tid.tsv"),read(p.labels)+"excluded_999_b\tgenome_x2\ttid_999_excluded.fa\t/assemblies/x2.fa\t"+repeat('e')+"\tcfg\tserver\tref1\tunknown\td__unknown;p__unknown\tunknown\tunknown\ttraining_exclusion\t"+repeat('5')+"\n");
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.build(p.counts.toString(),p.roster.toString(),p.layout.toString(),duplicateExcludedLabel.toString(),p.exclusions.toString(),dir.resolve("bad-label-tid.tsv").toString());}},"duplicate excluded label tid");
		String unknownDomain="unknown\td__unknown;p__unknown\tunknown\tunknown\tno_match";
		String badUnknown=read(p.labels).replace(uc+"\tgenome_"+uc+"\ttid_103_c.fa\t/assemblies/c.fa\t"+repeat('c')+"\tcfg\tserver\tref1\t"+unknownDomain,uc+"\tgenome_"+uc+"\ttid_103_c.fa\t/assemblies/c.fa\t"+repeat('c')+"\tcfg\tserver\tref1\tunknown\td__unknown;p__unknown\tBacteria\tunknown\tno_match");
		Path badUnknownPath=write(dir.resolve("labels-unknown-domain.tsv"),badUnknown);
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.build(p.counts.toString(),p.roster.toString(),p.layout.toString(),badUnknownPath.toString(),p.exclusions.toString(),dir.resolve("bad-unknown-domain.tsv").toString());}},"unknown label domain");
		String badClassified=read(p.labels).replace("\tBacteria\tFirmicutes\t-\t","\tArchaea\tFirmicutes\t-\t");
		Path badClassifiedPath=write(dir.resolve("labels-classified-domain.tsv"),badClassified);
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.build(p.counts.toString(),p.roster.toString(),p.layout.toString(),badClassifiedPath.toString(),p.exclusions.toString(),dir.resolve("bad-classified-domain.tsv").toString());}},"classified domain disagreement");
		String badPartial=read(p.labels).replace("\tArchaea\tunknown\tno_phylum\t","\tBacteria\tunknown\tno_phylum\t");
		Path badPartialPath=write(dir.resolve("labels-partial-domain.tsv"),badPartial);
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.build(p.counts.toString(),p.roster.toString(),p.layout.toString(),badPartialPath.toString(),p.exclusions.toString(),dir.resolve("bad-partial-domain.tsv").toString());}},"partial domain disagreement");
		Path neitherRostered=write(dir.resolve("labels-neither-rostered.tsv"),read(p.labels)+"g_bogus\tgenome_bogus\ttid_1000_bogus.fa\t/assemblies/bogus.fa\t"+repeat('e')+"\tcfg\tserver\tref1\tunknown\td__unknown;p__unknown\tunknown\tunknown\tfixture\t"+repeat('6')+"\n");
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.build(p.counts.toString(),p.roster.toString(),p.layout.toString(),neitherRostered.toString(),p.exclusions.toString(),dir.resolve("bad-neither-rostered.tsv").toString());}},"label neither rostered nor excluded");
		Path malformedMetadata=write(dir.resolve("table-malformed-metadata.tsv"),read(out).replace("#item_count\t7","#item_count\tnot-a-count"));
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.load(malformedMetadata.toString(),p.counts.toString(),p.roster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString());}},"malformed persisted table metadata");
		String[] lines=read(out).split("\n",-1);int first=-1,second=-1;for(int i=0;i<lines.length;i++){if(lines[i].length()>0&&!lines[i].startsWith("#")){if(first<0){first=i;}else{second=i;break;}}}check(first>=0&&second>=0,"persisted table has two payload rows");String swap=lines[first];lines[first]=lines[second];lines[second]=swap;final Path malformedOrder=write(dir.resolve("table-malformed-item-order.tsv"),join(lines));
		expectFailure(new Runnable(){public void run(){MagQCExpectedCopyTableTyped.load(malformedOrder.toString(),p.counts.toString(),p.roster.toString(),p.layout.toString(),p.labels.toString(),p.exclusions.toString());}},"malformed persisted table item order");
	}

	private static Paths inputs(Path dir) throws Exception{
		Path family=write(dir.resolve("families.tsv"),"#rank\trep_id\tocc_total\n0\tf0\t3\n1\tf1\t2\n");
		String familyHash=sha(family);
		String ra="tid_101_a.fa",rb="tid_102_b.fa",rc="tid_103_c.fa",ua=unit(ra),ub=unit(rb),uc=unit(rc);
		check("g_b205e3669b467ab2903a1bec2239a348e3f45204b24055c92812610c73d943a7".equals(ua),"manifest LF unit digest");
		Path layout=write(dir.resolve("layout.tsv"),"#schema_version\ttyped_item_layout_v1\n#family_list\t"+family+"\t#\tsha256\t"+familyHash+"\n#columns\t"+MagQCExpectedCopyTableTyped.LAYOUT_COLUMNS+"\n"
				+"P\t0\tf0\t0\nP\t1\tf1\t1\nN\tr16\t-\t2\nN\tr23\t-\t3\nN\tr5\t-\t4\nN\trother\t-\t5\nN\ttrna\t-\t6\n");
		Path exclusions=write(dir.resolve("exclusions.tsv"),"#excluded_tid\tclass\treason\n999\ttraining\tfixture\n");
		String rosterHeader="#schema_version\tclean_assembly_roster_v1\n#columns\t"+MagQCExpectedCopyTableTyped.ROSTER_COLUMNS+"\n";
		Path roster=write(dir.resolve("roster.tsv"),rosterHeader
				+ua+"\tgenome_"+ua+"\t"+ra+"\t/assemblies/a.fa\t"+repeat('a')+"\t101\tBacteria\n"
				+ub+"\tgenome_"+ub+"\t"+rb+"\t/assemblies/b.fa\t"+repeat('b')+"\t102\tArchaea\n"
				+uc+"\tgenome_"+uc+"\t"+rc+"\t/assemblies/c.fa\t"+repeat('c')+"\t103\tBacteria\n");
		String labelHeader="#columns\t"+MagQCExpectedCopyTableTyped.LABEL_COLUMNS+"\n";
		Path labels=write(dir.resolve("labels.tsv"),labelHeader
				+ua+"\tgenome_"+ua+"\t"+ra+"\t/assemblies/a.fa\t"+repeat('a')+"\tcfg\tserver\tref1\tclassified\td__Bacteria;p__Firmicutes\tBacteria\tFirmicutes\t-\t"+repeat('1')+"\n"
				+ub+"\tgenome_"+ub+"\t"+rb+"\t/assemblies/b.fa\t"+repeat('b')+"\tcfg\tserver\tref1\tpartial\td__Archaea;p__unknown\tArchaea\tunknown\tno_phylum\t"+repeat('2')+"\n"
				+uc+"\tgenome_"+uc+"\t"+rc+"\t/assemblies/c.fa\t"+repeat('c')+"\tcfg\tserver\tref1\tunknown\td__unknown;p__unknown\tunknown\tunknown\tno_match\t"+repeat('3')+"\n"
				+"excluded_999\tgenome_x\ttid_999_excluded.fa\t/assemblies/x.fa\t"+repeat('d')+"\tcfg\tserver\tref1\tunknown\td__unknown;p__unknown\tunknown\tunknown\ttraining_exclusion\t"+repeat('4')+"\n");
		Path counts=write(dir.resolve("counts.tsv"),"#schema_version\tclean_typed_counts_v1\n#columns\t"+MagQCExpectedCopyTableTyped.COUNT_COLUMNS+"\n"
				+ub+"\tgenome_"+ub+"\t102\tP\t0\t4\n"
				+ua+"\tgenome_"+ua+"\t101\tP\t0\t2\n"
				+ua+"\tgenome_"+ua+"\t101\tP\t1\t1\n"
				+ub+"\tgenome_"+ub+"\t102\tN\tr16\t2\n"
				+uc+"\tgenome_"+uc+"\t103\tN\tr16\t0\n"
				+ua+"\tgenome_"+ua+"\t101\tN\tr16\t1\n"
				+ua+"\tgenome_"+ua+"\t101\tN\tr23\t1\n"
				+ub+"\tgenome_"+ub+"\t102\tN\tr23\t0\n"
				+uc+"\tgenome_"+uc+"\t103\tN\tr23\t0\n"
				+ua+"\tgenome_"+ua+"\t101\tN\tr5\t0\n"
				+ub+"\tgenome_"+ub+"\t102\tN\tr5\t0\n"
				+uc+"\tgenome_"+uc+"\t103\tN\tr5\t0\n"
				+ua+"\tgenome_"+ua+"\t101\tN\trother\t0\n"
				+ub+"\tgenome_"+ub+"\t102\tN\trother\t0\n"
				+uc+"\tgenome_"+uc+"\t103\tN\trother\t0\n"
				+ua+"\tgenome_"+ua+"\t101\tN\ttrna\t1\n"
				+ub+"\tgenome_"+ub+"\t102\tN\ttrna\t0\n"
				+uc+"\tgenome_"+uc+"\t103\tN\ttrna\t0\n");
		return new Paths(family,layout,roster,labels,exclusions,counts);
	}

	private static final class Paths{final Path family,layout,roster,labels,exclusions,counts;Paths(Path f,Path l,Path r,Path q,Path x,Path c){family=f;layout=l;roster=r;labels=q;exclusions=x;counts=c;}}
	private static String read(Path p) throws Exception{return new String(Files.readAllBytes(p),StandardCharsets.UTF_8);}
	private static Path write(Path p,String s) throws Exception{Files.write(p,s.getBytes(StandardCharsets.UTF_8));return p;}
	private static String unit(String sourceRel) throws Exception{byte[] d=MessageDigest.getInstance("SHA-256").digest((sourceRel+"\n").getBytes(StandardCharsets.UTF_8));StringBuilder b=new StringBuilder(66);b.append("g_");for(byte x:d){b.append(String.format("%02x",x&255));}return b.toString();}
	private static String repeat(char c){StringBuilder b=new StringBuilder(64);for(int i=0;i<64;i++){b.append(c);}return b.toString();}
	private static String sha(Path p) throws Exception{byte[] d=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(p));StringBuilder b=new StringBuilder(64);for(byte x:d){b.append(String.format("%02x",x&255));}return b.toString();}
	private static void expectFailure(Runnable r,String what){try{r.run();throw new RuntimeException("FAIL: did not reject "+what);}catch(IllegalArgumentException expected){}}
	private static void expectFailureContaining(Runnable r,String what,String fragment){try{r.run();throw new RuntimeException("FAIL: did not reject "+what);}catch(IllegalArgumentException expected){check(expected.getMessage()!=null&&expected.getMessage().indexOf(fragment)>=0,"specific failure for "+what+": "+expected.getMessage());}}
	private static String join(String[] lines){StringBuilder b=new StringBuilder();for(int i=0;i<lines.length;i++){if(i>0){b.append('\n');}b.append(lines[i]);}return b.toString();}
	private static boolean close(double a,double b){return Math.abs(a-b)<=1e-12;}
	private static void check(boolean ok,String message){if(!ok){throw new RuntimeException("FAIL: "+message);}}
	private static void deleteTree(Path p) throws Exception{if(!Files.exists(p)){return;}java.nio.file.DirectoryStream<Path> ds=Files.newDirectoryStream(p);try{for(Path child:ds){if(Files.isDirectory(child)){deleteTree(child);}else{Files.deleteIfExists(child);}}}finally{ds.close();}Files.deleteIfExists(p);}
}
