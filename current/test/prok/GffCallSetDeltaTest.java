package test.prok;

import prok.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.charset.StandardCharsets;

/** Multiplicity, family attribution, score changes and sortedness fixtures.
 * @author Raiden, Ganyu
 */
public final class GffCallSetDeltaTest {
	public static void main(String[] args) throws Exception{
		final Path root=Files.createTempDirectory("gff_delta_");final String a="a\tBBTools\tCDS\t1\t120\t50\t+\t0\tCDS\n",b="b\tBBTools\trRNA\t1\t120\t50\t+\t0\t5S,model:euk5S_animal,family:euk5S\n",c=b.replace("\t50\t","\t51\t"),d="c\tBBTools\trRNA\t1\t120\t50\t+\t0\t5S\n";
		// Actual caller emits legacy attributes followed by commas.
		final String legacy=d.replace("\t5S\n","\t5S,len:120\n");
		final String its="d\tBBTools\tinternal_transcribed_spacer_region\t121\t220\t.\t+\t0\tproduct:ITS,from:a,to:b,len:100\n";
		final String its1=its.replace("d\tBBTools","e\tBBTools").replace("product:ITS,","product:ITS1,");
		final String its2=its.replace("d\tBBTools","f\tBBTools").replace("product:ITS,","product:ITS2,");
		Files.write(root.resolve("a.gff"),(a+a+b+legacy+its+its1+its2).getBytes(StandardCharsets.US_ASCII));Files.write(root.resolve("b.gff"),(a+c+legacy+its+its1+its2.replace("\t220\t","\t221\t")).getBytes(StandardCharsets.US_ASCII));final GffCallSetDelta tool=new GffCallSetDelta();tool.run(root.resolve("a.gff").toString(),root.resolve("b.gff").toString(),root.resolve("result").toString());
		check(tool.counts.get("CDS").removed==1 && tool.counts.get("CDS").shared==1,"Duplicate emitted calls must retain multiset multiplicity");check(tool.counts.get("euk5S").removed==1 && tool.counts.get("euk5S").added==1,"Changed score must remain an exact-row change within the same family");check(tool.counts.get("legacy5S").shared==1,"Legacy5S must remain distinct from euk5S");
		check(tool.counts.get("ITS").shared==1 && tool.counts.get("ITS1").shared==1,"Derived spacer products must be conserved separately");
		check(tool.counts.get("ITS2").removed==1 && tool.counts.get("ITS2").added==1,"Changed spacer boundaries must remain visible as full-row differences");
		Files.write(root.resolve("unsorted.gff"),(b+a).getBytes(StandardCharsets.US_ASCII));boolean rejected=false;try{new GffCallSetDelta().run(root.resolve("unsorted.gff").toString(),root.resolve("b.gff").toString(),root.resolve("bad").toString());}catch(IllegalArgumentException e){rejected=e.getMessage().contains("byte-sorted");}check(rejected,"Incorrect sort order must fail rather than produce false deltas");
		System.out.println("PASS GffCallSetDeltaTest multiplicity family score_change sortedness");
		bindings(root);
	}
	static void bindings(Path root)throws Exception{
		final String a="a\tBBTools\tCDS\t1\t120\t50\t+\t0\tCDS";
		final String b=a.replace("a\tBBTools","b\tBBTools"),missing=a.replace("a\tBBTools","missing\tBBTools"),u=a.replace("a\tBBTools","u\tBBTools");
		final String rows="REMOVED\t"+a+"\nADDED\t"+a+"\nADDED\t"+b+"\nADDED\t"+missing+"\nREMOVED\t"+a+"\nREMOVED\t"+u+'\n';
		final Path changes=root.resolve("binding.changes.tsv"),loci=root.resolve("binding.loci.tsv"),ungraded=root.resolve("binding.ungraded.tsv"),audit=root.resolve("binding.audit.tsv");
		Files.write(changes,(GffCallSetDelta.CHANGE_HEADER+'\n'+rows).getBytes(StandardCharsets.UTF_8));
		Files.write(loci,(GffCallSetDelta.PANEL_LOCI_HEADER+'\n'+owner("before","A","a")+owner("after","B","a")
			+owner("after","C","b")+owner("after","D","b")+owner("after","C","b")+owner("SOURCE","ignored","a")).getBytes(StandardCharsets.UTF_8));
		Files.write(ungraded,(GffCallSetDelta.PANEL_UNGRADED_HEADER+"\nbefore\tF\tPROK\t1\tu\tBBTools\tRNA\t1\t120\t50\t+\t0\tRNA\n").getBytes(StandardCharsets.UTF_8));
		final String prefix=root.resolve("bound").toString();
		GffCallSetDelta.bind(changes.toString(),loci.toString(),ungraded.toString(),null,null,prefix);
		final String result=read(prefix+".bindings.tsv.gz"),summary=read(prefix+".binding_summary.tsv");
		check(summary.contains("UNIQUE\t4\n") && summary.contains("AMBIGUOUS\t1\n") && summary.contains("UNRESOLVED\t1\n"),
			"All six native delta rows, including duplicates, need explicit ownership outcomes");
		check(result.contains("UNIQUE\t1\tA\tPANEL_ANNOTATION_NATIVE_ROWS\tREMOVED\t"+a)
			&& result.contains("UNIQUE\t1\tB\tPANEL_ANNOTATION_NATIVE_ROWS\tADDED\t"+a),"Ownership must stay arm-specific");
		check(result.contains("AMBIGUOUS\t2\tC,D\t") && result.contains("UNIQUE\t1\tF\t"),"Repeated metadata rows cannot inflate owners; ungraded rows supply ownership too");
		final StringBuilder original=new StringBuilder();for(String line:result.split("\n")){if(line.startsWith("binding_status")){continue;}final String[] fields=line.split("\t",5);check(fields.length==5,"Binding adds exactly four fields");original.append(fields[4]).append('\n');}
		check(original.toString().equals(rows),"Binding must retain every byte of each original delta row, in original order");
		Files.write(audit,(NcrnaEndpointAudit.HEADER+'\n'+auditOwner("before","P1","A","a")+auditOwner("after","P1","B","a")
			+auditOwner("after","P2","wrong_panel","a")+auditOwner("after","P1","C","b")+auditOwner("after","P1","D","b")
			+auditOwner("before","P1","F","u")).getBytes(StandardCharsets.UTF_8));
		GffCallSetDelta.bind(changes.toString(),null,null,audit.toString(),"P1",root.resolve("audit_bound").toString());
		final String bound=read(root.resolve("audit_bound.bindings.tsv.gz").toString());
		check(!bound.contains("wrong_panel") && bound.contains("UNIQUE\t1\tB\tFROZEN_AUDIT_P1"),"Explicit panel prevents a shared seqid in another cohort from becoming an owner");
		Files.write(root.resolve("bad_owners.tsv"),"wrong\theader\n".getBytes(StandardCharsets.UTF_8));
		boolean rejected=false;try{GffCallSetDelta.bind(changes.toString(),root.resolve("bad_owners.tsv").toString(),ungraded.toString(),null,null,root.resolve("bad_binding").toString());}catch(IllegalArgumentException e){rejected=true;}
		check(rejected,"Unknown metadata schema must fail loudly");
		System.out.println("GFF_DELTA_BINDING_TEST_PASS arm_panel_isolation unique_ambiguous_unresolved owners_deduplicated raw_row_multiplicity_preserved schema_rejection");
	}
	static String owner(String arm,String genome,String seq){return arm+'\t'+genome+"\tCDS\tCALL\t"+seq+"\t1\t120\t+\t0\t1\tTP\t1\t1\t120\t0\t0\t0\tLINEAR\n";}
	static String auditOwner(String arm,String panel,String genome,String seq){return arm+'\t'+panel+'\t'+genome+'\t'+seq+"\tPROKARYOTIC\tSRP\tCALL\tcall\tSRP\t1\t120\t+\tTP\ttruth\t0\t0\tFULLSPAN_INTACT\tFIXTURE\n";}
	static String read(String path){final fileIO.ByteFile in=fileIO.ByteFile.makeByteFile(path,false);final structures.ByteBuilder b=new structures.ByteBuilder();
		try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{check(!in.close(),"Read full fixture output");}return new String(b.toBytes(),StandardCharsets.UTF_8);}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
