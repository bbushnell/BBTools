package test.prok;

import java.io.File;
import java.util.Arrays;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.*;
import prok.NcrnaClippingComparison.Genome;
import prok.NcrnaClippingComparison.Contig;
import prok.NcrnaClippingComparison.Feature;
import structures.ByteBuilder;

/** Actual native-schema pipeline plus assembly-isolation, duplicates and6S subtype/union controls.
 * @author Brian, Raiden */
public final class NcrnaFreshPanelGradeTest {
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="),"Use outdir=fresh_directory");
		final String d=args[0].substring(7);check(new File(d).mkdir(),"Fresh fixture required");
		final String e=d+"/evidence",n=e+"/native",proposal=d+"/proposal",loci=d+"/loci",candidate=d+"/candidate",disposition=d+"/disposition";
		ProkPanelCmEvidenceTest.main(new String[]{"outdir="+e});
		EukPanelCmProposal.main(new String[]{"prokinputs=t","inputs="+n+"/good.assemblies.tsv","sources="+n+"/sources.tsv","contigs="+n+"/good.contigs.tsv","contigsources="+n+"/sources.tsv","verified="+e+"/verified.tsv",
			"hits="+e+"/pass.hits.tsv.gz","pairs="+e+"/pass.hit_overlaps.tsv.gz","summary="+e+"/pass.summary.tsv","features="+e+"/pass.sources.tsv.gz","associations="+e+"/pass.associations.tsv.gz",
			"expected=2","expectedmodels=2","expectedruns=2","expectedhits=3","expectedsources=7","expectedassociations=3","outprefix="+proposal});
		EukPanelLoci.main(new String[]{"in="+proposal+".overlay.tsv.gz","models="+e+"/models.tsv","scans="+proposal+".scan_denominators.tsv","expectedhits=3","expectedmodels=2","expectedscans=4","outprefix="+loci});
		check(rows(loci+".invalid_hits.tsv.gz")==1 && rows(loci+".locus_members.tsv.gz")==2,"One invalid hit and two valid members conserve all three raw hits");
		EukPanelLocusDiagnostics.main(new String[]{"loci="+loci+".loci.tsv.gz","members="+loci+".locus_members.tsv.gz","expectedmembers=2","outprefix="+loci});
		EukPanelSourceMask.main(new String[]{"in="+proposal+".source_review.tsv.gz","sources="+n+"/sources.tsv","expectedrows=4","outprefix="+loci});
		EukPanelTruthCandidate.main(new String[]{"loci="+loci+".loci.tsv.gz","members="+loci+".locus_members.tsv.gz","links="+proposal+".source_links.tsv.gz","sources="+n+"/sources.tsv","models="+e+"/models.tsv",
			"cases="+loci+".locus_cases.tsv","diagnostics="+loci+".diagnostic_cases.tsv","expectedloci=2","expectedmodels=2","expectedmembers=2","expectedlinks=3","expectedambiguous=0","expectedcrossfold=0","expectedunexpected=0","outprefix="+candidate});
		EukPanelTruthDisposition.main(new String[]{"policy="+POLICY,"candidate="+candidate+".truth_candidate.tsv.gz","members="+loci+".locus_members.tsv.gz","cases="+loci+".locus_cases.tsv","expectedloci=2","expectedmembers=2","expectedambiguous=0","outprefix="+disposition});
		check(text(loci+".locus_members.tsv.gz").equals(text(disposition+".locus_members.tsv.gz")),"Fresh prok disposition must not absorb or rewrite members");
		EukPanelTruthFreeze.main(new String[]{"policy="+POLICY,"in="+disposition+".truth_candidate.tsv.gz","out="+d+"/truth.tsv.gz","expectedloci=2"});
		check(text(d+"/truth.tsv.gz").equals(text(disposition+".truth_candidate.tsv.gz").replace("FREEZE_CANDIDATE_NOT_FROZEN","TRUTH_V1")),"Fresh freeze changes status only");
		boolean rejected=false;try{EukPanelTruthFreeze.main(new String[]{"in="+disposition+".truth_candidate.tsv.gz","out="+d+"/wrong_policy.tsv","expectedloci=2"});}catch(IllegalArgumentException ex){rejected=true;}check(rejected,"Legacy default must not silently accept fresh-prok policy");
		EukPanelHistoricalMask.main(new String[]{"truth="+d+"/truth.tsv.gz","review="+proposal+".source_review.tsv.gz","links="+proposal+".source_links.tsv.gz","models="+e+"/models.tsv","sources="+n+"/sources.tsv","outprefix="+d+"/history","expectedloci=2","expectedmodels=2","expectedreview=4","expectedlinks=3"});
		check(rows(d+"/history.historical_masks.tsv.gz")==5,"All five historical RNA rows survive, including invalid/unresolved rows");
		final StringBuilder verified=new StringBuilder(text(e+"/verified.tsv")),runs=new StringBuilder(NcrnaFreshPanelGrade.RUN_HEADER+'\n');
		for(String file:new String[]{d+"/truth.tsv.gz",d+"/history.historical_masks.tsv.gz",loci+".out_of_map.tsv.gz",d+"/runs.tsv"}){verified.append(HASH).append('\t').append(file).append('\n');}
		for(int i=1;i<=2;i++){
			final String r=d+"/caller"+i;check(new File(r).mkdir(),"Fresh saved caller fixture");
			write(r+"/calls.gff",i==1?call("shared",1,8,"+","prok5S")+call("shared",1,8,"+","prok5S")+call("tid|11|tagged",1,4,"-","tRNA"):call("shared",1,4,"+","prok5S"));
			write(r+"/INPUTS.sha80",HASH+'\t'+n+(i==1?"/a.fa\n":"/b.fa\n"));
			write(r+"/time.txt","wall_seconds=1 user_seconds=1 system_seconds=0\n");
			write(r+"/stderr","ncRNA Alignments (rnasep): 0\n"+(i==1?"tRNA Out: 1\n":"tRNA Alignments: 0\n")+"ncRNA Shared Sweep Passes: 2\nncRNA Aligned Cells: 0\n");
			write(r+"/MANIFEST.sha80","Synthetic reviewed saved-run manifest pointer\n");
			for(String file:new String[]{"calls.gff","INPUTS.sha80","time.txt","stderr","MANIFEST.sha80"}){verified.append(HASH).append('\t').append(r).append('/').append(file).append('\n');}
			runs.append(i).append("\tGCF_").append(i).append(".1\t").append(r).append('\t').append(HASH).append('\t').append(r).append("/calls.gff\t").append(HASH).append('\n');
		}
		write(d+"/verified.tsv",verified.toString());write(d+"/runs.tsv",runs.toString());run(d,e,n,loci,d+"/grade");
		final String work=text(d+"/grade.workload.tsv");
		check(work.contains("SHIPPED\tbacteria\t1\talignments\ttRNA\tNA\n") && work.contains("SHIPPED\tarchaea\t2\talignments\ttRNA\t0\n"),"Unprinted tRNA alignment work is NA; a measured zero remains zero");
		workloadChecks(d);
		final ByteFile summary=RrnaResourceIO.open(d+"/grade.summary.tsv");final LineParser1 p=new LineParser1('\t');int checked=0;
		try{for(byte[] line;(line=summary.nextLine())!=null;){if(line[0]=='#'){continue;}p.set(line);if(!p.termEquals("5S",2)){continue;}
			if(p.termEquals("bacteria",1)){check(p.parseLong(3)==1 && p.parseLong(4)==1 && p.parseLong(9)==1,"Duplicate is FP, not second TP");checked++;}
			if(p.termEquals("archaea",1)){check(p.parseLong(3)==0 && p.parseLong(4)==1,"Shared seqid must never match another assembly's truth");checked++;}
			if(p.termEquals("all",1)){check(p.parseLong(3)==1 && p.parseLong(4)==2 && p.parseLong(5)==0,"Pooled counters must conserve domain outcomes");checked++;}
		}}finally{check(!summary.close(),"Read summary");}check(checked==3,"Both domain and pooled5S rows required");
		write(d+"/runs.tsv",runs.toString().replace("2\tGCF_2.1","2\tGCF_1.1"));rejected=false;
		try{run(d,e,n,loci,d+"/swapped");}catch(IllegalArgumentException ex){check(ex.getMessage().contains("selected index"),"Wrong swapped-run failure");rejected=true;}check(rejected,"Wrong assembly cannot borrow another run");
		sixs(d);spanning(d);srp(d);batchGff(d);
		System.out.println("FRESH_PROK_TRUTH_GRADE_TEST_PASS native_schema_conservation no_euk_absorption explicit_policy assembly_isolation duplicates sixs_subtypes_union partial_masks endpoints srp_physical_union subtype_agreement");
	}
	static void batchGff(String d){
		final String root=d+"/batch_path";check(new File(root).mkdir(), "Fresh batch-path fixture required");
		final String plain=root+"/calls.gff", compressed=root+"/calls.gff.gz";
		final String rows=call("chr", 1, 8, "+", "prok5S")+call("chr", 1, 8, "+", "prok5S")+call("chr", 20, 30, "-", "tRNA");
		write(plain, rows);write(compressed, rows);
		check(NcrnaFreshPanelGrade.batchCallsPath(root, null).equals(plain)
			&& NcrnaFreshPanelGrade.batchCallsPath(root, compressed).equals(compressed), "Default raw GFF path is unchanged and explicit compressed suffix is retained");
		boolean rejected=false;try{NcrnaFreshPanelGrade.batchCallsPath(root, d+"/other/calls.gff.gz");}catch(IllegalArgumentException ex){rejected=true;}
		check(rejected, "A batch cannot borrow another run's GFF");
		rejected=false;try{NcrnaFreshPanelGrade.batchCallsPath(null, compressed);}catch(IllegalArgumentException ex){rejected=true;}
		check(rejected, "An override without batchroot must fail");
		java.util.ArrayList<Feature> before=null;
		for(int i=0;i<2;i++){
			final NcrnaClippingComparison grade=new NcrnaClippingComparison();final Genome g=new Genome("A", "bacteria", "/unused", HASH);
			final Contig c=new Contig(g, "chr", 100, "PROKARYOTIC");grade.addGenome(g);g.contigs.put(c.id, c);g.ordered.add(c);
			final ByteStreamWriter unknown=RrnaResourceIO.writer(root+"/unknown"+i+".tsv");
			try{grade.loadCalls(NcrnaFreshPanelGrade.batchCallsPath(root, i==0 ? null : compressed), "SHIPPED", g.contigs, unknown, new ByteBuilder());}
			finally{RrnaResourceIO.finish(unknown);}
			check(c.calls.size()==3, "Both compression paths retain every row, including duplicate calls");
			if(i==0){before=new java.util.ArrayList<Feature>(c.calls);}
			else{for(int j=0;j<before.size();j++){final Feature a=before.get(j), z=c.calls.get(j);
				final NcrnaGffFamilyGrader.Locus al=a.locus(a.family), zl=z.locus(z.family);
				check(al.start==zl.start && al.stop==zl.stop && al.strand==zl.strand && al.family.equals(zl.family) && al.seqid.equals(zl.seqid), "Compressed and raw streams supply identical ordered features to grading");}}
		}
		System.out.println("FRESH_BATCH_GFF_PATH_TEST_PASS default_raw explicit_compressed source_root_binding native_reader_parity duplicates");
	}
	static void run(String d,String e,String n,String loci,String out){
		NcrnaFreshPanelGrade.main(new String[]{"inputs="+n+"/good.assemblies.tsv","sources="+n+"/sources.tsv","contigs="+n+"/good.contigs.tsv","contigsources="+n+"/sources.tsv","verified="+d+"/verified.tsv","models="+e+"/models.tsv",
			"truth="+d+"/truth.tsv.gz","history="+d+"/history.historical_masks.tsv.gz","masks="+loci+".out_of_map.tsv.gz","runs="+d+"/runs.tsv","expected=2","expectedmodels=2","expectedloci=2","expectedhistoryrows=5","outprefix="+out});
	}
	static void workloadChecks(String d){
		final ByteStreamWriter out=RrnaResourceIO.writer(d+"/workload_negative.tsv");
		try{
			boolean rejected=false;try{NcrnaClippingComparison.readWorkload(d+"/caller1/stderr",out,new ByteBuilder(),"TEST","bacteria","1");}catch(IllegalArgumentException ex){rejected=true;}
			check(rejected,"Legacy workload contract still requires its tRNA alignment counter");
			write(d+"/incomplete.stderr","ncRNA Alignments (rnasep): 1\nncRNA Shared Sweep Passes: 2\n");
			rejected=false;try{NcrnaClippingComparison.readWorkload(d+"/incomplete.stderr",out,new ByteBuilder(),"TEST","bacteria","1",false);}catch(IllegalArgumentException ex){rejected=true;}
			check(rejected,"Optional tRNA accounting cannot hide a missing global aligned-cell counter");
		}finally{RrnaResourceIO.finish(out);}
	}
	static void sixs(String d){
		final String a=sixsTruth(1,100,"RF00013",true,"FULLSPAN_INTACT"),z=sixsTruth(2,300,"RF01685",true,"PARTIAL_INTERNAL"),mask=sixsTruth(3,500,"RF00013",false,"PARTIAL_INTERNAL");
		write(d+"/sixs.tsv",EukPanelTruthFreeze.HEADER+'\n'+a+z+mask);write(d+"/empty_history.tsv",NcrnaClippingComparison.HISTORY_HEADER+'\n');
		final NcrnaClippingComparison x=new NcrnaClippingComparison();final Genome g=new Genome("A","bacteria","/unused",HASH);x.addGenome(g);
		final Contig c=new Contig(g,"chr",1000,"PROKARYOTIC");g.contigs.put(c.id,c);g.ordered.add(c);x.modelFamilies.put("RF00013","6S");x.modelFamilies.put("RF01685","6S");
		x.loadFrozenTruth(d+"/sixs.tsv",3,POLICY);x.loadHistoricalMasks(d+"/empty_history.tsv",0);
		c.calls.add(new Feature(c,"6S","sixs_rf00013",100,199,"+","q1","NA","CALL"));c.calls.add(new Feature(c,"6S","sixs_rf00013",100,199,"+","duplicate","NA","CALL"));
		c.calls.add(new Feature(c,"6S","sixs_rf00013",300,399,"+","wrong_subtype","NA","CALL"));c.calls.add(new Feature(c,"6S","sixs_rf01685",500,599,"+","other_mask_subtype","NA","CALL"));
		c.masks.add(new Feature(c,"6S","6S",700,799,"+","unresolved_history","NA","HIST_ONLY"));
		c.calls.add(new Feature(c,"6S","sixs_rf00013",700,799,"+","unknown_history_subtype","NA","CALL"));
		final ByteStreamWriter audit=RrnaResourceIO.writer(d+"/sixs.audit.tsv");try{for(String view:new String[]{"6S_RF00013","6S_RF01685","6S"}){x.grade(c,"TEST",view,audit,new ByteBuilder());}}finally{RrnaResourceIO.finish(audit);}
		check(Arrays.equals(x.totals.get("TEST\tbacteria\tPROKARYOTIC\t6S_RF00013"),new long[]{1,3,0,1,1,1,1,0,1,1}),"Wrong-subtype/duplicate calls remain FP; generic6S history cannot supply subtype-specific evidence");
		check(Arrays.equals(x.totals.get("TEST\tbacteria\tPROKARYOTIC\t6S_RF01685"),new long[]{0,1,1,0,0,0,0,0,0,0}),"Wrong-subtype call misses true subtype and cannot use another subtype's below-GA mask");
		check(Arrays.equals(x.totals.get("TEST\tbacteria\tPROKARYOTIC\t6S"),new long[]{2,1,0,2,2,2,1,2,2,2}),"Union retains physical matches and both below-GA and unresolved-history masks without inventing subtype support");
		check(x.truthStrata.get("TEST\tbacteria\tPROKARYOTIC\t6S\tPARTIAL_INTERNAL")[0]==1,"Partial internal truths remain reported, not silently filtered");
	}
	static String sixsTruth(int id,int start,String accession,boolean admitted,String geometry){
		final String[] p=NcrnaFrozenLocusComparisonTest.truth(id,start,start+99,admitted,false,geometry).trim().split("\t",-1);
		p[6]="6S";p[7]="bacteria";p[8]="PROKARYOTIC";p[9]=accession;p[19]=accession;p[25]=admitted?accession:"NA";p[34]=POLICY;return String.join("\t",p)+'\n';
	}
	static void spanning(String d){
		final int[] starts={100,100,300},ends={399,199,399},bits={300,100,90};final String[] models={"RF00177","RF01959","RF00177"};
		final StringBuilder candidate=new StringBuilder(EukPanelTruthDisposition.INPUT_HEADER+'\n'),members=new StringBuilder(EukPanelTruthDisposition.MEMBER_HEADER+'\n');
		for(int i=0;i<3;i++){
			candidate.append(EukPanelTruthDispositionTest.locus(i+1,starts[i],ends[i],"SSU",models[i],bits[i],i==0?"SPANNING_AMBIGUOUS":"NONE",false).replace("\ttest\tNUCLEAR\t","\tbacteria\tPROKARYOTIC\t"));
			members.append(EukPanelTruthDispositionTest.member(i+1,starts[i],ends[i],"SSU",models[i],bits[i]).replace("\tNUCLEAR\t","\tPROKARYOTIC\t").replace("\ttest\t1000\t","\tbacteria\t1000\t"));
		}
		write(d+"/spanning_candidate.tsv",candidate.toString());write(d+"/spanning_members.tsv",members.toString());
		write(d+"/spanning_cases.tsv",EukPanelTruthDisposition.CASE_HEADER+"\n1\tSPANNING_SEED\tA\tchr\t+\tSSU\trun\t1\t100\t399\trun\t2\t100\t199\trun\t3\t300\t399\n");
		EukPanelTruthDisposition.main(new String[]{"policy="+POLICY,"candidate="+d+"/spanning_candidate.tsv","members="+d+"/spanning_members.tsv","cases="+d+"/spanning_cases.tsv","expectedloci=3","expectedmembers=3","expectedambiguous=1","outprefix="+d+"/spanning"});
		check(rows(d+"/spanning.truth_candidate.tsv.gz")==3 && rows(d+"/spanning.merge_map.tsv")==0,"Flagged prok spanning parent cannot absorb either witness locus");
		check(text(d+"/spanning.locus_members.tsv.gz").equals(members.toString()),"Every original prok member stays under its original locus");
		final String[] old=candidate.toString().split("\n"),now=text(d+"/spanning.truth_candidate.tsv.gz").split("\n");
		for(int i=1;i<old.length;i++){check(now[i].startsWith(old[i]+'\t'),"Fresh disposition preserves every original candidate field and uncertainty flag");}
	}
	static void srp(String d){
		final String prefix=d+"/srp_original",out=d+"/srp_union";
		final String[] raw={EukPanelLociTest.row(1,"RF01854","SRP_large",100,399,200,true),EukPanelLociTest.row(2,"RF00169","SRP_small",200,299,100,true),
			EukPanelLociTest.row(3,"RF00169","SRP_small",600,699,60,true),EukPanelLociTest.row(4,"RF00005","tRNA",850,949,70,true)};
		final StringBuilder hits=new StringBuilder(EukPanelLoci.HEADER+'\n');for(String row:raw){hits.append(row.replace("\ttest\t","\tbacteria\t").replace("\tUNKNOWN\t","\tPROKARYOTIC\t")).append('\n');}
		write(d+"/srp_hits.tsv",hits.toString());
		write(d+"/srp_models.tsv","bundle_sha80\taccession\tfamily\tsubtype\tclen\tga\tsource_cm\tsource_cm_sha80\n"+EukPanelLociTest.model("RF01854","SRP_large")+EukPanelLociTest.model("RF00169","SRP_small")+EukPanelLociTest.model("RF00005","tRNA"));
		write(d+"/srp_scans.tsv","selection_bin\t"+EukPanelCmProposal.SUMMARY_HEADER+'\n'+(EukPanelLociTest.scan("RF01854",1)+EukPanelLociTest.scan("RF00169",2)+EukPanelLociTest.scan("RF00005",1)).replace("test\tA","bacteria\tA"));
		EukPanelLoci.main(new String[]{"in="+d+"/srp_hits.tsv","models="+d+"/srp_models.tsv","scans="+d+"/srp_scans.tsv","outprefix="+prefix,"expectedhits=4","expectedmodels=3","expectedscans=3"});
		NcrnaSrpUnionTruth.main(new String[]{"loci="+prefix+".loci.tsv.gz","members="+prefix+".locus_members.tsv.gz","expectedmembers=4","expectedoldsrp=3","outprefix="+out});
		check(rows(out+".truth.tsv.gz")==2 && rows(out+".members.tsv.gz")==3,"Nested small/large hits become one physical locus while the separate small locus survives");
		check(text(out+".summary.tsv").contains("3\t2\t3\t4\t2\t0\t0\t"+NcrnaSrpUnionTruth.POLICY),"Union member and old/new locus counts conserve the declared SRP subset");
		final String originalMembers=text(prefix+".locus_members.tsv.gz");for(String member:text(out+".members.tsv.gz").split("\n")){if(member.startsWith("assembly\t")){continue;}
			check(originalMembers.contains(member.substring(0,member.lastIndexOf('\t'))+'\n'),"Every union member retains its complete original43-field row");}
		final StringBuilder frozen=new StringBuilder(EukPanelTruthFreeze.HEADER+'\n');final NcrnaClippingComparison x=new NcrnaClippingComparison();
		final Genome g=new Genome("A","bacteria","/unused",HASH);x.addGenome(g);final Contig c=new Contig(g,"chr",1000,"PROKARYOTIC");g.contigs.put(c.id,c);g.ordered.add(c);
		for(String line:text(prefix+".loci.tsv.gz").split("\n")){if(line.startsWith("locus_id\t")){continue;}final String[] fields=line.split("\t");
			final String[] t=NcrnaFrozenLocusComparisonTest.truth(Integer.parseInt(fields[0]),Integer.parseInt(fields[4]),Integer.parseInt(fields[5]),true,false,"FULLSPAN_INTACT").trim().split("\t",-1);
			t[6]=fields[6];t[7]="bacteria";t[8]="PROKARYOTIC";t[9]=fields[9];t[19]=fields[9];t[25]=fields[9];t[34]=POLICY;frozen.append(String.join("\t",t)).append('\n');x.modelFamilies.put(fields[9],fields[6]);}
		write(d+"/srp_frozen_v1.tsv",frozen.toString());x.loadFrozenTruth(d+"/srp_frozen_v1.tsv",4,POLICY);
		Feature control=null;for(Feature f:c.truth){if(f.family.equals("tRNA")){control=f;}}
		NcrnaSrpUnionTruth.apply(x,out+".truth.tsv.gz",out+".old_ids.tsv",3,2);
		check(c.truth.size()==3 && c.truth.contains(control),"Non-SRP truth object is retained unchanged, only3old SRP rows replaced by2physical rows");
		c.calls.add(new Feature(c,"SRP_small","srp_small",100,399,"+","wrong_label","NA","CALL"));
		c.calls.add(new Feature(c,"SRP_small","srp_small",600,699,"+","right_label","NA","CALL"));
		final ByteStreamWriter audit=RrnaResourceIO.writer(d+"/srp_audit.tsv"),agreement=RrnaResourceIO.writer(d+"/srp_agreement.tsv");
		try{for(String view:new String[]{"SRP_small","SRP_large","SRP"}){x.grade(c,"TEST",view,audit,new ByteBuilder());}
			check(Arrays.equals(NcrnaFreshPanelGrade.srpAgreement(c,agreement,new ByteBuilder()),new long[]{2,1}),"Physical detection and subtype agreement must remain separate");
		}finally{RrnaResourceIO.finish(audit,agreement);}
		final long[] small=x.totals.get("TEST\tbacteria\tPROKARYOTIC\tSRP_small"),large=x.totals.get("TEST\tbacteria\tPROKARYOTIC\tSRP_large"),union=x.totals.get("TEST\tbacteria\tPROKARYOTIC\tSRP");
		check(small[0]==1 && small[1]==1 && small[2]==0 && large[0]==0 && large[2]==1 && union[0]==2 && union[1]==0 && union[2]==0,"Wrong subtype can be a physical union TP while remaining subtype FP/FN");
		boolean rejected=false;final LineParser1 p=new LineParser1('\t');final byte[] invalid=raw[3].getBytes(java.nio.charset.StandardCharsets.UTF_8);p.set(invalid);
		try{new EukPanelLoci.Hit(invalid,p,"SRP_UNION");}catch(IllegalArgumentException ex){rejected=true;}check(rejected,"Union override cannot regroup a different biological family");
		final EukPanelLoci.Hit full=unionHit(EukPanelLociTest.row(10,"RF01854","SRP_large",100,399,90,true)),partial=unionHit(EukPanelLociTest.partial(EukPanelLociTest.row(11,"RF00169","SRP_small",250,299,100,true)));
		final java.util.ArrayList<EukPanelLoci.Hit> group=new java.util.ArrayList<EukPanelLoci.Hit>();group.add(partial);group.add(full);
		final java.util.ArrayList<EukPanelLoci.Locus> physical=new java.util.ArrayList<EukPanelLoci.Locus>();EukPanelLoci.buildPartition(group,physical);
		check(physical.size()==1 && physical.get(0).seed==full && EukPanelLoci.subtype(physical.get(0).members)==partial,"Physical seed class and highest-bit subtype stay independent across the two original SRP folds");
	}
	static EukPanelLoci.Hit unionHit(String row){final byte[] raw=row.getBytes(java.nio.charset.StandardCharsets.UTF_8);final LineParser1 p=new LineParser1('\t');p.set(raw);return new EukPanelLoci.Hit(raw,p,"SRP_UNION");}
	static String call(String seq,int a,int z,String strand,String family){return seq+"\tBBTools\t"+(family.equals("tRNA")?"tRNA":"rRNA")+'\t'+a+'\t'+z+"\t.\t"+strand+"\t.\tRNA,family:"+family+'\n';}
	static long rows(String path){return text(path).split("\n").length-1;}
	static String text(String path){return EukPanelCmProposalTest.read(path);}
	static void write(String path,String value){
		// Negative fixtures deliberately replace their own inputs; product outputs stay fresh.
		final ByteStreamWriter out=new ByteStreamWriter(path,true,false,false);out.start();out.print(value);RrnaResourceIO.finish(out);
	}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	static final String HASH="00000000000000000000",POLICY=EukPanelTruthDisposition.PROK_POLICY;
}
