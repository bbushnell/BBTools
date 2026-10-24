package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMitoCmEvidence;
import prok.NcrnaMitoCmEvidence.Hit;
import prok.NcrnaMitoCmEvidence.Span;
import parse.LineParser1;
import structures.ByteBuilder;

/** Coordinate examples use explicit hand-derived genomic bases, not min/max shortcuts.
 * @author Ganyu, Raiden */
public final class NcrnaMitoCmEvidenceTest {
	public static void main(String[] args){
		require(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require outdir=NEW_DIRECTORY sha80=CANONICAL_EXECUTABLE");
		final String root=args[0].substring(7);require(new File(root).mkdir(),"Fixture output must be new");
		final String sha80=args[1].substring(6);require(new File(sha80).canExecute(),"Canonical sha80 is required for real fixture input bindings");
		final ArrayList<Span> join=spans("91..100:+,1..10:+");
		check(NcrnaMitoCmEvidence.project(join,8,14,'+'),"98..100:+,1..4:+","Forward joined hit must cross the origin through both segments");
		check(NcrnaMitoCmEvidence.project(join,14,8,'-'),"1..4:-,98..100:-","Opposite candidate hit reverses segment order and genomic strands");
		final ArrayList<Span> complement=spans("1..10:-,91..100:-");
		check(NcrnaMitoCmEvidence.project(complement,3,14,'+'),"1..8:-,97..100:-","Complement join must project through the actual transcript order");
		check(NcrnaMitoCmEvidence.project(complement,14,3,'-'),"97..100:+,1..8:+","Reverse hit on complement join returns plus genomic segments");
		final ArrayList<Span> gapped=NcrnaMitoCmEvidence.project(spans("11..20:+,41..50:+"),8,13,'+');
		check(gapped,"18..20:+,41..43:+","Unaligned genomic gap cannot become hit evidence");
		require(NcrnaMitoCmEvidence.overlap(NcrnaMitoCmEvidence.normalize(gapped),spans("21..40:+"),true)==0,"Min/max projection would incorrectly cover the twenty-base join gap");
		final ArrayList<Span> circular=NcrnaMitoCmEvidence.genomicSpans(98,104,'+',100,10);
		check(circular,"98..100:+,1..4:+","Circular genomic window preserves exact origin-spanning footprint");
		check(NcrnaMitoCmEvidence.genomicSpans(104,98,'-',100,10),"1..4:-,98..100:-","Circular reverse genomic window keeps orientation");
		final ArrayList<Span> original=NcrnaMitoCmEvidence.genomicSpans(2,5,'+',100,10),appended=NcrnaMitoCmEvidence.genomicSpans(102,105,'+',100,10);
		final String key=NcrnaMitoCmEvidence.physicalKey("genomic","RF00177","A",original);
		require(key.equals(NcrnaMitoCmEvidence.physicalKey("genomic","RF00177","A",appended)),"Appended and original representations share one exact physical key");
		require(!key.equals(NcrnaMitoCmEvidence.physicalKey("candidate","RF00177","A",appended)),"Target universes must remain distinct");
		require(!key.equals(NcrnaMitoCmEvidence.physicalKey("genomic","RF01960","A",appended)),"Different model scores cannot collapse into one winner");
		check(appended,"2..5:+","A lone valid appended hit retains a physical footprint rather than disappearing");
		final Hit copy=new Hit(),origin=new Hit();copy.physicalGroupKey=key;copy.appended=true;copy.score=100;copy.mf=1;copy.mt=4;copy.line=3;
		origin.physicalGroupKey=key;origin.score=90;origin.mf=1;origin.mt=4;origin.line=2;
		require(NcrnaMitoCmEvidence.chooseRepresentative(null,copy)==copy,"Lone appended representation is retained");
		require(NcrnaMitoCmEvidence.chooseRepresentative(copy,origin)==origin && NcrnaMitoCmEvidence.chooseRepresentative(origin,copy)==origin,"Exact duplicate representations prefer original coordinates independent of input order");
		final Hit otherUniverse=new Hit();otherUniverse.physicalGroupKey=NcrnaMitoCmEvidence.physicalKey("candidate","RF00177","A",appended);otherUniverse.score=150;
		badGeometry(()->NcrnaMitoCmEvidence.chooseRepresentative(origin,otherUniverse));
		final ArrayList<Span> united=NcrnaMitoCmEvidence.normalize(spans("1..5:+,3..8:+"));
		check(united,"1..8:+","Overlapping source segments form a union for physical overlap counts");
		require(NcrnaMitoCmEvidence.overlap(united,spans("4..7:+"),true)==4 && NcrnaMitoCmEvidence.overlap(united,spans("4..7:-"),false)==4,"Same and opposite strand overlaps count physical bases once");
		badGeometry(()->NcrnaMitoCmEvidence.project(join,1,21,'+'));
		badGeometry(()->NcrnaMitoCmEvidence.project(join,14,8,'+'));
		badGeometry(()->NcrnaMitoCmEvidence.genomicSpans(1,105,'+',100,10));
		badGeometry(()->NcrnaMitoCmEvidence.parseSegments("1..101:+",100));
		final String row="A - bacterial_SSU RF00177 cm 1 10 9 18 + no 1 0.50 0.0 65.0 1e-8 ! original description";
		final String footer="# Program: cmsearch\n# Version: 1.1.5 (Sep 2023)\n# [ok]\n";
		write(root+"/hit.tbl","#target name accession query name\n"+row+"\n"+footer);
		final ArrayList<Hit> hits=new ArrayList<Hit>();NcrnaMitoCmEvidence.parseTable(root+"/hit.tbl","RF00177",100,50,hits::add);
		require(hits.size()==1 && hits.get(0).from==9 && hits.get(0).to==18 && hits.get(0).score==65 && hits.get(0).evalue.equals("1e-8"),"Native parser must preserve the exact raw coordinates and scientific-notation E-value");
		write(root+"/zero.tbl","#target name accession query name\n"+footer);
		final int[] count={0};NcrnaMitoCmEvidence.parseTable(root+"/zero.tbl","RF00177",100,50,h->{count[0]++;});require(count[0]==0,"A complete zero-hit scan remains a valid run");
		final String[] bad={row+"\n"+footer.replace("# [ok]\n",""),row+"\n"+footer.replace("1.1.5","1.0.0"),
			row.replace("RF00177","RF01960")+"\n"+footer,row.replace(" cm "," hmm ")+"\n"+footer,
			row.replace("9 18 +","18 9 +")+"\n"+footer,row.replace("cm 1 10","cm 0 10")+"\n"+footer,
			row.replace("65.0","NaN")+"\n"+footer,row.replace(" + no "," + unknown ")+"\n"+footer,
			row+"\n"+footer+row+"\n","A - too few\n"+footer,row.replace("1e-8","-1e-8")+"\n"+footer};
		for(int i=0;i<bad.length;i++){final String path=root+"/bad_"+i+".tbl";write(path,bad[i]);boolean failed=false;
			try{NcrnaMitoCmEvidence.parseTable(path,"RF00177",100,50,h->{});}catch(IllegalArgumentException e){failed=true;}
			require(failed,"Malformed/incomplete table must fail, case "+i);
		}
		badGeometry(()->NcrnaMitoCmEvidence.execution(null,null,"ln001.jgi"));
		badGeometry(()->NcrnaMitoCmEvidence.execution("123","executor=head pid=456","ln001.jgi"));
		badGeometry(()->NcrnaMitoCmEvidence.execution(null,"executor=head pid=0","ln001.jgi"));
		badGeometry(()->NcrnaMitoCmEvidence.execution("NA",null,"ln001.jgi"));
		badGeometry(()->NcrnaMitoCmEvidence.execution(null,"executor=head pid=456",""));
		endToEnd(root,sha80,false);
		endToEnd(root,sha80,true);
		for(String suffix:new String[]{".raw_rows.tsv.gz",".loci.tsv",".genomes.tsv",".universal_pairs.tsv",".cross.tsv",".cost.tsv"}){
			require(read(root+"/cli/ok"+suffix).equals(read(root+"/cli_head/ok"+suffix)),"Executor metadata must not change scientific output: "+suffix);
		}
		System.out.println("MITO_CM_EXECUTION_TEST_PASS slurm=true head=true scientific_parity=true invalid_metadata=5 missing_binding=true");
		partitioned(root,sha80);
		restrictedGenomic(root,sha80);
		System.out.println("MITO_CM_EVIDENCE_TEST_PASS joined_complement=true circular_arbitration=true universe_separation=true zero_hit=true cross_cost=true cli_raw=4 cli_groups=3 cli_runs=4 cli_negatives=2 geometry_negatives=5 table_negatives=11");
	}
	private static void partitioned(String root,String sha80){
		final Fixture f=new Fixture(root+"/cli_partitioned",sha80);require(new File(f.root).mkdir(),"Partition fixture must be new");f.parts=2;f.membership="membership.tsv";
		for(String name:new String[]{"selection.tsv","loci.tsv","members.tsv","pairs.tsv"}){f.put(name,read(root+"/cli/"+name));}
		final String header="universe\tpartition\taccession\ttarget\toriginal_length\tsearch_length\n";
		final String members=header+"genomic\t0\tA\tA\t100\t110\ngenomic\t1\tB\tB\t100\t100\ncandidate\t0\tA\tA_rrna_1\t20\t20\n";
		f.put("membership.tsv",members);
		final ByteBuilder models=new ByteBuilder("accession\tpath\tsha80\n");int job=1000;
		for(String model:new String[]{"RF00177","RF01960"}){
			f.put(model+".cm",read(root+"/cli/"+model+".cm"));models.append(model).tab().append(f.root+"/"+model+".cm").tab().append(f.hash(f.root+"/"+model+".cm")).nl();
			for(String kind:new String[]{"genomic","candidate"}){job++;for(int part=0;part<2;part++){
				final String dir="runs/"+model+"_"+kind+"_p"+part;require(new File(f.root+"/"+dir).mkdirs(),"Partition run must be new");
				final String original=root+"/cli/runs/"+model+"_"+kind;
				f.put(dir+"/hits.tbl",part==0?read(original+"/hits.tbl"):"#target name\n# Program: cmsearch\n# Version: 1.1.5 (Sep 2023)\n# [ok]\n");
				f.put(dir+"/job_id.txt",job+"\n");
				final String[] files={"hits.tbl","job_id.txt","command.txt","time.txt","host.txt","INPUTS.sha80","stdout","stderr","producer.sh"};
				for(int i=2;i<files.length;i++){f.put(dir+"/"+files[i],read(original+"/"+files[i]));}
				final ByteBuilder manifest=new ByteBuilder();for(String name:files){manifest.append(f.hash(f.root+"/"+dir+"/"+name)).tab().append(name).nl();}
				f.put(dir+"/MANIFEST.sha80",manifest.toString());f.put(dir+"/COMPLETE.sha80",f.hash(f.root+"/"+dir+"/MANIFEST.sha80")+"\n");f.put(dir+"/exit_code.txt","0\n");
			}}
		}
		f.put("models.tsv",models.toString());f.verified("verified.sha80",null);f.run("ok","loci.tsv","verified.sha80");
		for(String suffix:new String[]{".loci.tsv",".genomes.tsv",".universal_pairs.tsv",".cross.tsv"}){require(read(root+"/cli/ok"+suffix).equals(read(f.root+"/ok"+suffix)),"Partitioning must preserve scientific summaries: "+suffix);}
		require(lines(read(f.root+"/ok.runs.tsv"))==9 && lines(read(f.root+"/ok.raw_rows.tsv.gz"))==5,"All eight partition runs, including six zero-hit runs, and all four raw hits must survive");
		final String execution=read(f.root+"/ok.execution.tsv");require(execution.contains("RF00177_genomic_p0\tslurm\t1001\tln001.jgi\tNA\n") && execution.contains("RF00177_genomic_p1\tslurm\t1001\tln001.jgi\tNA\n"),"Two tasks retain their actual shared scheduler allocation without inventing job IDs");
		final ByteFile in=ByteFile.makeByteFile(f.root+"/ok.cost.tsv",true);final LineParser1 p=new LineParser1('\t');
		try{require(new String(in.nextLine(),StandardCharsets.UTF_8).contains("sum_task_wall_seconds"),"Packed allocation cost must name task wall, not fleet elapsed or allocation wall");for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.parseInt(3)==4 && p.parseDouble(4)==4 && p.parseDouble(8)==4,"All four task times per universe must count");require(p.parseInt(2)==(p.termEquals("genomic",0)?210:20),"Partition cost must count full target bases only once");}}finally{require(!in.close(),"Partition cost read failed");}
		final String[] bad={members.replace("genomic\t0\tA", "genomic\t1\tA"),members+"candidate\t0\tA\tA_rrna_1\t20\t20\n",header+"genomic\t0\tA\tA\t100\t110\n",members.replace("genomic\t1\tB","genomic\t2\tB")};
		for(int i=0;i<bad.length;i++){f.membership="bad_membership_"+i+".tsv";f.put(f.membership,bad[i]);f.verified("bad_verified_"+i+".sha80",null);f.reject("bad_membership_"+i,"loci.tsv","bad_verified_"+i+".sha80");}
		f.membership="membership.tsv";f.verified("missing_part.sha80",f.root+"/runs/RF01960_candidate_p1/hits.tbl");f.reject("missing_part","loci.tsv","missing_part.sha80");
		f.verified("unbound_membership.sha80",f.root+"/membership.tsv");f.reject("unbound_membership","loci.tsv","unbound_membership.sha80");
		System.out.println("MITO_CM_PARTITION_TEST_PASS runs=8 zero_hit_runs=6 raw_hits=4 scientific_parity=true shared_scheduler_id=true bases_once=true negatives=6");
	}
	/** A has two usable annotations and only candidate searches; unannotated B has a zero-hit genomic search. */
	private static void restrictedGenomic(String root,String sha80){
		final Fixture f=new Fixture(root+"/cli_restricted",sha80);require(new File(f.root).mkdir(),"Restricted fixture must be new");
		f.parts=2;f.lociCount=2;f.membership="membership.tsv";f.genomicParts=1;f.genomicBases=100;f.eligibility="eligibility.tsv";f.genomicMembership="genomic_membership.tsv";
		f.put("selection.tsv",read(root+"/cli/selection.tsv"));
		f.put("loci.tsv",read(root+"/cli/loci.tsv")+"A_rrna_2\tA\t1\tmitochondrion\t2\t41..60\t5S\tEXTRACTED\tfalse\tfalse\tNA\t+\t41\t60\t41..60:+\tfalse\t20\tfixture.gbff\n");
		f.put("members.tsv",read(root+"/cli/members.tsv")+"A_rrna_2\tseq_0\t20\n");f.put("pairs.tsv",read(root+"/cli/pairs.tsv"));
		final String header="universe\tpartition\taccession\ttarget\toriginal_length\tsearch_length\n";
		final String membership=header+"genomic\t0\tA\tA\t100\t110\ngenomic\t1\tB\tB\t100\t100\ncandidate\t0\tA\tA_rrna_1\t20\t20\ncandidate\t0\tA\tA_rrna_2\t20\t20\n";
		f.put("membership.tsv",membership);f.put("genomic_membership.tsv",header+"genomic\t0\tB\tB\t100\t100\n");
		final String eligibility="accession\trrna_features\tusable_rrna\tgenomic_selected\toriginal_length\tprefix_added\tsearch_length\nA\t2\t2\tfalse\t100\t10\t110\nB\t0\t0\ttrue\t100\t0\t100\n";
		f.put("eligibility.tsv",eligibility);final ByteBuilder models=new ByteBuilder("accession\tpath\tsha80\n");int job=2000;
		for(String model:new String[]{"RF00177","RF01960"}){
			f.put(model+".cm",read(root+"/cli/"+model+".cm"));models.append(model).tab().append(f.root+'/'+model+".cm").tab().append(f.hash(f.root+'/'+model+".cm")).nl();
			copyFixtureRun(f,"genomic/"+model+"_genomic_p0",root+"/cli_partitioned/runs/"+model+"_genomic_p1",++job);
			for(int part=0;part<2;part++){copyFixtureRun(f,"runs/"+model+"_candidate_p"+part,root+"/cli_partitioned/runs/"+model+"_candidate_p"+part,++job);}
		}
		f.put("models.tsv",models.toString());f.verified("verified.sha80",null);f.run("ok","loci.tsv","verified.sha80");
		require(lines(read(f.root+"/ok.runs.tsv"))==7,"Two models times one genomic plus two candidate partitions require six completed runs");
		final String loci=read(f.root+"/ok.loci.tsv"),genomes=read(f.root+"/ok.genomes.tsv"),cross=read(f.root+"/ok.cross.tsv");
		require(loci.contains("A_rrna_1\tA\tseq_0\tNA\tNA\t1\t0\t0\tSSU_SUPPORT_REVIEW") && loci.contains("\tNOT_SEARCHED\n"),"An unsearched genome must have NA genomic evidence while preserving its candidate support");
		require(genomes.contains("\tNOT_SEARCHED\tSEARCHED\n") && genomes.contains("\tSEARCHED\tNO_ANNOTATED_TARGET\n"),"Unsearched genomes and searched zero-hit genomes must have distinct coverage labels");
		require(cross.contains("NOT_SEARCHED\tSSU_SUPPORT_REVIEW\tTIE\t1\t1\t1\t0\t2\n") && cross.contains("NOT_SEARCHED\tNO_SAME_STRAND_GA_HIT_UNRESOLVED\tTIE\t1\t1\t1\t0\t2\n"),"Evidence/preference cross-table must stratify genomic coverage and conserve both candidate annotations");
		final String cost=read(f.root+"/ok.cost.tsv");require(cost.contains("genomic\t1\t100\t2\t") && cost.contains("candidate\t2\t40\t4\t"),"Cost denominators reflect each actual searched universe, not the full census for both");
		f.put("bad_eligibility.tsv",eligibility.replace("A\t2\t2\tfalse","A\t2\t2\ttrue"));f.eligibility="bad_eligibility.tsv";f.verified("bad_eligibility.sha80",null);f.reject("bad_eligibility","loci.tsv","bad_eligibility.sha80");f.eligibility="eligibility.tsv";
		f.put("missing_genomic.tsv",header);f.genomicMembership="missing_genomic.tsv";f.verified("missing_genomic.sha80",null);f.reject("missing_genomic","loci.tsv","missing_genomic.sha80");
		f.put("excluded_genomic.tsv",header+"genomic\t0\tB\tB\t100\t100\ngenomic\t0\tA\tA\t100\t110\n");f.genomicMembership="excluded_genomic.tsv";f.verified("excluded_genomic.sha80",null);f.reject("excluded_genomic","loci.tsv","excluded_genomic.sha80");f.genomicMembership="genomic_membership.tsv";
		f.put("missing_candidate.tsv",membership.replace("candidate\t0\tA\tA_rrna_2\t20\t20\n",""));f.membership="missing_candidate.tsv";f.verified("missing_candidate.sha80",null);f.reject("missing_candidate","loci.tsv","missing_candidate.sha80");f.membership="membership.tsv";
		f.genomicBases=101;f.reject("wrong_genomic_bases","loci.tsv","verified.sha80");
		System.out.println("MITO_CM_RESTRICTED_TEST_PASS source_eligibility=true separate_roots=true unequal_partitions=true unsearched_NA=true coverage_strata=true bases_once=true negatives=5");
	}
	private static void copyFixtureRun(Fixture f,String dir,String source,int job){
		require(new File(f.root+'/'+dir).mkdirs(),"Restricted fixture run must be new");
		final String[] files={"hits.tbl","command.txt","time.txt","host.txt","INPUTS.sha80","stdout","stderr","producer.sh","job_id.txt"};
		for(String name:files){f.put(dir+'/'+name,name.equals("job_id.txt")?job+"\n":read(source+'/'+name));}
		final ByteBuilder manifest=new ByteBuilder();for(String name:files){manifest.append(f.hash(f.root+'/'+dir+'/'+name)).tab().append(name).nl();}
		f.put(dir+"/MANIFEST.sha80",manifest.toString());f.put(dir+"/COMPLETE.sha80",f.hash(f.root+'/'+dir+"/MANIFEST.sha80")+"\n");f.put(dir+"/exit_code.txt","0\n");
	}
	private static void endToEnd(String root,String sha80,boolean head){
		final Fixture f=new Fixture(root+(head?"/cli_head":"/cli"),sha80);require(new File(f.root).mkdir(),"CLI fixture directory must be new");
		f.put("selection.tsv","accession\ttaxid\torganelle_class\tannotation_group\tlength_quartile0\toriginal_length\tprefix_added\tsearch_length\tpartial\tpseudo\tcontig_fallback\tcanary\n"
			+"A\t1\tmitochondrion\tannotated\t0\t100\t10\t110\tfalse\tfalse\tfalse\ttrue\nB\t2\tmitochondrion\tno_rrna\t0\t100\t0\t100\tfalse\tfalse\tfalse\tfalse\n");
		final String locusHeader="locus_id\taccession\ttaxid\torganelle_class\trrna_index\traw_location\tproduct\tstatus\tpartial\tpseudo\tpseudogene\tstrand\tmin_coordinate\tmax_coordinate\tsegments_in_transcript_order\tcrosses_origin\textracted_length\tgbff\n";
		final String locus=locusHeader+"A_rrna_1\tA\t1\tmitochondrion\t1\tcomplement(join(91..100,1..10))\tSSU\tEXTRACTED\tfalse\tfalse\tNA\t-\t1\t100\t1..10:-,91..100:-\ttrue\t20\tfixture.gbff\n";
		f.put("loci.tsv",locus);
		f.put("members.tsv","locus_id\tunique_id\tlength\nA_rrna_1\tseq_0\t20\n");
		f.put("pairs.tsv","unique_id\tmodel\tcandidate_length\tmodel_length\tquery_side\tquantum_identity\ttrace_identity\tref_start0\tref_stop0\tmatches\tsubstitutions\tinsertions\tdeletions\tneutral_N\tpaired_columns\tcandidate_paired_fraction\tmodel_paired_fraction\n"
			+"seq_0\tuniversal\t20\t20\tcandidate\t1\t1\t0\t19\t20\t0\t0\t0\t0\t20\t1\t1\nseq_0\teuk18S_universal\t20\t20\tcandidate\t1\t1\t0\t19\t20\t0\t0\t0\t0\t20\t1\t1\n");
		final String[] models={"RF00177","RF01960"};final ByteBuilder modelRows=new ByteBuilder("accession\tpath\tsha80\n");int job=1000;
		for(String model:models){
			// Metadata-only unit fixture, never passed to cmsearch or presented as a real CM.
			f.put(model+".cm","INFERNAL1/a [synthetic metadata fixture]\nACC "+model+"\nCLEN 20\nGA 50\n");
			modelRows.append(model).tab().append(f.root+"/"+model+".cm").tab().append(f.hash(f.root+"/"+model+".cm")).nl();
			for(String kind:new String[]{"genomic","candidate"}){
				final String dir="runs/"+model+"_"+kind;require(new File(f.root+"/"+dir).mkdirs(),"Synthetic run directory must be new");
				String rows="";
				if(model.equals("RF00177")){rows=kind.equals("genomic")?
					"A - SSU RF00177 cm 1 10 104 95 - no 1 0.5 0 70 1e-8 ! crossing\nA - SSU RF00177 cm 1 6 8 3 - no 1 0.5 0 90 1e-10 ! original\nA - SSU RF00177 cm 1 6 108 103 - no 1 0.5 0 100 1e-12 ! appended\n":
					"A_rrna_1 - SSU RF00177 cm 1 12 3 14 + no 1 0.5 0 80 1e-9 ! joined\n";}
				f.put(dir+"/hits.tbl","#target name\n"+rows+"# Program: cmsearch\n# Version: 1.1.5 (Sep 2023)\n# [ok]\n");
				f.put(dir+"/command.txt","cmsearch --cpu 1 --cut_ga --noali synthetic_fixture\n");
				f.put(dir+"/time.txt","wall_seconds=1 user_seconds=1 system_seconds=0\n");
				final String executionFile=head?"execution.txt":"job_id.txt";
				f.put(dir+"/"+executionFile,(head?"executor=head pid=":"")+(++job)+"\n");
				for(String name:new String[]{"INPUTS.sha80","stdout","stderr","producer.sh","CODE.txt","SCOPE.txt"}){f.put(dir+"/"+name,"synthetic unit fixture\n");}
				f.put(dir+"/host.txt","ln001.jgi\n");
				final ByteBuilder manifest=new ByteBuilder();
				for(String name:new String[]{"hits.tbl","command.txt","time.txt",executionFile,"INPUTS.sha80","stdout","stderr","producer.sh","CODE.txt","host.txt","SCOPE.txt"}){manifest.append(f.hash(f.root+"/"+dir+"/"+name)).tab().append(name).nl();}
				f.put(dir+"/MANIFEST.sha80",manifest.toString());f.put(dir+"/COMPLETE.sha80",f.hash(f.root+"/"+dir+"/MANIFEST.sha80")+"\n");f.put(dir+"/exit_code.txt","0\n");
			}
		}
		f.put("models.tsv",modelRows.toString());f.verified("verified.sha80",null);
		f.run("ok","loci.tsv","verified.sha80");
		final String r=f.root+"/ok";
		require(read(r+".execution.tsv").contains(head?"RF00177_genomic\thead\tNA\tln001.jgi\t1001\n":"RF00177_genomic\tslurm\t1001\tln001.jgi\tNA\n"),"Head PID and scheduler ID must occupy separate explicit fields");
		require(lines(read(r+".raw_rows.tsv.gz"))==5 && lines(read(r+".hits.tsv.gz"))==5,"Full CLI must preserve all four raw hits");
		require(lines(read(r+".runs.tsv"))==5 && lines(read(r+".genomes.tsv"))==3 && lines(read(r+".loci.tsv"))==2 && lines(read(r+".universal_pairs.tsv"))==3,"Full CLI conserves four runs including zero-hit arms, two genomes, one annotation and both universal rows");
		require(read(r+".genomes.tsv").contains("B\t2\tmitochondrion\tno_rrna\t0\t0\tNO_GA_HIT_UNRESOLVED\n"),"No-annotation/no-hit genome must remain explicitly unresolved");
		require(read(r+".loci.tsv").contains("A_rrna_1\tA\tseq_0\t2\t0\t1\t0\t0\tSSU_SUPPORT_REVIEW\t"),"Exact duplicate collapse must yield two genomic and one candidate same-strand support groups");
		final String associations=read(r+".associations.tsv.gz");
		require(lines(associations)==5 && associations.contains("\tA_rrna_1\t10\t0\t20\t10\n") && associations.contains("\tA_rrna_1\t12\t0\t20\t12\n"),"Circular/complement/join associations must use exact10/12-base footprints");
		final ByteFile hitsIn=ByteFile.makeByteFile(r+".hits.tsv.gz",true);final LineParser1 p=new LineParser1('\t');int representatives=0,appendedCopies=0;
		try{hitsIn.nextLine();for(byte[] line;(line=hitsIn.nextLine())!=null;){p.set(line);if(p.termEquals("true",21)){representatives++;}if(p.parseInt(13)==108){require(p.termEquals("true",19) && p.termEquals("false",21),"Appended duplicate stays in raw output but is not the representative");appendedCopies++;}}}finally{require(!hitsIn.close(),"CLI hit read failed");}
		require(representatives==3 && appendedCopies==1,"Four retained rows must yield exactly three representatives");
		final ByteFile runsIn=ByteFile.makeByteFile(r+".runs.tsv",true);int zeroRuns=0;
		try{runsIn.nextLine();for(byte[] line;(line=runsIn.nextLine())!=null;){p.set(line);if(p.parseInt(9)==0){require(p.parseInt(10)==0,"Zero raw hits cannot create a physical group");zeroRuns++;}}}finally{require(!runsIn.close(),"CLI run read failed");}
		require(zeroRuns==2,"Both zero-hit model arms must be present");
		require(read(r+".cross.tsv").contains("SSU_SUPPORT_REVIEW\tTIE\t1\t1\t1\t0\t1\n"),"Cross-table must retain the sole SSU-supported annotation as an exact universal tie");
		final ByteFile costIn=ByteFile.makeByteFile(r+".cost.tsv",true);int costRows=0;
		try{costIn.nextLine();for(byte[] line;(line=costIn.nextLine())!=null;){p.set(line);require(p.parseInt(3)==2 && p.parseDouble(4)==2 && p.parseDouble(5)==1 && p.parseDouble(8)==2,"Two synthetic one-second jobs sum to2 wall/CPU but have maxwall1, not parallel elapsed2");
			if(p.termEquals("genomic",0)){require(p.parseInt(1)==2 && p.parseInt(2)==210,"Genomic cost denominator includes original100+100 plus10 circular prefix bases");}
			else{require(p.termEquals("candidate",0) && p.parseInt(1)==1 && p.parseInt(2)==20,"Candidate cost denominator is the extracted20-base transcript");}costRows++;}}
		finally{require(!costIn.close(),"CLI cost read failed");}require(costRows==2,"Both target universes need separate cost rows");
		f.put("missing_join.tsv",locus.replace("1..10:-,91..100:-","NA"));f.verified("bad_join_verified.sha80",null);
		f.reject("bad_join","missing_join.tsv","bad_join_verified.sha80");
		f.verified("missing_run_verified.sha80",f.root+"/runs/RF01960_candidate/hits.tbl");f.reject("missing_run","loci.tsv","missing_run_verified.sha80");
		f.verified("missing_execution_verified.sha80",f.root+"/runs/RF00177_genomic/"+(head?"execution.txt":"job_id.txt"));f.reject("missing_execution","loci.tsv","missing_execution_verified.sha80");
	}
	private static final class Fixture {
		Fixture(String r,String s){root=r;sha80=s;}
		void put(String name,String contents){final String path=root+"/"+name;write(path,contents);require(!paths.contains(path),"Fixture source path written twice");paths.add(path);}
		String hash(String path){
			final File output=new File(root+"/.hash.stdout");
			try{final Process process=new ProcessBuilder(sha80,path).redirectOutput(output).redirectError(ProcessBuilder.Redirect.INHERIT).start();require(process.waitFor()==0,"Canonical sha80 failed for fixture input");}
			catch(Exception e){throw new IllegalStateException("Canonical sha80 fixture binding failed",e);}
			final String value=read(output.getPath()).trim();require(value.matches("[0-9a-f]{20}"),"Canonical sha80 must return exactly twenty hex characters");return value;
		}
		void verified(String name,String excluded){final ByteBuilder b=new ByteBuilder();for(String path:paths){if(!path.equals(excluded)){b.append(hash(path)).tab().append(path).nl();}}write(root+"/"+name,b.toString());}
		void run(String name,String loci,String verified){final ArrayList<String> args=new ArrayList<String>(java.util.Arrays.asList("selection="+root+"/selection.tsv","loci="+root+"/"+loci,"models="+root+"/models.tsv","runsroot="+root+"/runs","members="+root+"/members.tsv","pairs="+root+"/pairs.tsv","verified="+root+"/"+verified,"outprefix="+root+"/"+name,"expectedgenomes=2","expectedloci="+lociCount,"expectedmodels=2","t=1","zl=9"));
			if(parts>0){args.add("partitions="+parts);args.add("membership="+root+"/"+membership);}
			if(eligibility!=null){args.add("genomicroot="+root+"/genomic");args.add("genomicpartitions="+genomicParts);args.add("genomicmembership="+root+'/'+genomicMembership);args.add("genomiceligibility="+root+'/'+eligibility);args.add("expectedgenomicgenomes=1");args.add("expectedgenomicbases="+genomicBases);}
			NcrnaMitoCmEvidence.main(args.toArray(new String[args.size()]));}
		void reject(String name,String loci,String verified){boolean failed=false;try{run(name,loci,verified);}catch(IllegalArgumentException e){failed=true;}require(failed && !new File(root+"/"+name+".runs.tsv").exists(),"Malformed join/missing declared run cannot produce a completed reduction");}
		final String root,sha80;int parts=0,lociCount=1,genomicParts=0,genomicBases=0;String membership,eligibility,genomicMembership;final ArrayList<String> paths=new ArrayList<String>();
	}
	private static String read(String path){final ByteFile in=ByteFile.makeByteFile(path,true);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{require(!in.close(),"Fixture read failed");}return new String(b.array,0,b.length(),StandardCharsets.UTF_8);}
	private static int lines(String s){int n=0;for(int i=0;i<s.length();i++){if(s.charAt(i)=='\n'){n++;}}return n;}
	private static ArrayList<Span> spans(String s){return NcrnaMitoCmEvidence.parseSegments(s,100);}
	private static void check(ArrayList<Span> spans,String expected,String why){require(NcrnaMitoCmEvidence.geometry(spans).equals(expected),why+": "+NcrnaMitoCmEvidence.geometry(spans));}
	private static void badGeometry(Runnable action){boolean failed=false;try{action.run();}catch(IllegalArgumentException e){failed=true;}require(failed,"Impossible geometry cannot be projected");}
	private static void write(String path,String text){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,true,false,false,false));w.start();w.print(text.getBytes(StandardCharsets.UTF_8));require(!w.poisonAndWait(),"Fixture write failed");}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
