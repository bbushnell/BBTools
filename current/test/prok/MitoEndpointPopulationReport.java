package test.prok;

import prok.RrnaResourceIO;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;
import tax.TaxNode;
import tax.TaxTree;

/** Frozen-taxonomy population report for already-scored validation anchors.
 * Operational algal roots are explicit inputs; unresolved hosts are retained.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden */
public final class MitoEndpointPopulationReport {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,MitoEndpointPopulationReport.class,false);
		try{final MitoEndpointPopulationReport x=new MitoEndpointPopulationReport();
			for(String arg:pp.args){final String[] kv=Parse.splitOnFirst(arg,'=');require(kv.length==2,"Expected explicit flag=value inputs");kv[0]=kv[0].toLowerCase(Locale.ROOT);
				if(kv[0].equals("evaluation")){x.evaluation=kv[1];}else if(kv[0].equals("enrollment")){x.enrollment=kv[1];}
				else if(kv[0].equals("gates")){x.gateFile=kv[1];}
				else if(kv[0].equals("tree")){x.treeFile=kv[1];}else if(kv[0].equals("rules")){x.ruleFile=kv[1];}
				else if(kv[0].equals("outprefix")){x.prefix=kv[1];}else{throw new IllegalArgumentException(kv[0]);}}
			require(x.evaluation!=null && x.enrollment!=null && x.treeFile!=null && x.ruleFile!=null && x.prefix!=null,"Require saved evaluation, enrollment, frozen taxonomy, explicit rules and fresh output");
			Shared.setThreads(1);x.run();
		}finally{Shared.closeStream(pp.outstream);}
	}
	private void run(){
		ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;ReadWrite.ZIPLEVEL=9;ReadWrite.setZipThreads(1);
		RrnaResourceIO.fresh(prefix,".members.tsv.gz",".summary.tsv",".offsets.tsv",".within.tsv",".roots.tsv",".onebase.tsv",".shift3plus1.tsv");
		if(gateFile!=null){RrnaResourceIO.fresh(prefix,".gate_summary.tsv",".gate_decisions.tsv",".gate_changes.tsv.gz");}
		final boolean old=TaxTree.LOAD_FAKE_TREE;try{TaxTree.LOAD_FAKE_TREE=false;tree=TaxTree.loadTaxTreeFile(treeFile);}finally{TaxTree.LOAD_FAKE_TREE=old;}
		loadRules();policyFixture();loadEnrollment();
		final ObjectMap<String,Stats> groups=new ObjectMap<String,Stats>(String.class,Stats.class);
		for(String name:new String[]{"ALL","CHORDATE","ARTHROPOD","FUNGAL","ALGAL_LINEAGE","OTHER"}){groups.put("GROUP\t"+name+"\tALL\tALL",new Stats());}
		final ByteFile in=RrnaResourceIO.open(evaluation);final LineParser1 p=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();
		final ByteStreamWriter detail=RrnaResourceIO.writer(prefix+".members.tsv.gz");int seen=0,scored=0;
		try{RrnaResourceIO.header(in,EVALUATION_HEADER,evaluation);
			detail.println("locus\thost_taxid\toriginal_phylum\tcurrent_phylum\tgroup\tmatched_ancestor_id\tmatched_ancestor_name\tassignment_reason\traw_error5\traw_error3\tnet_error5\tnet_error3\tlineage_leaf_to_root");
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==18,"Native evaluation schema changed");seen++;
				final Member m=members.get(p.parseString(0));require(m!=null && !m.seen && p.parseInt(1)==m.taxid,"Evaluation rows must bind each enrolled host/locus once");m.seen=true;
				if(p.termEquals("NOT_EVALUATED",17)){continue;}
				require(p.termEquals("VALIDATION",2) && p.termEquals("COMPLETE_SOURCE_ASSOCIATED",3) && p.parseInt(5)==192
					&& p.termEquals("VALIDATION_SEQUENTIAL_SCORER_NOT_WHOLE_GENOME_CALLER_ACCURACY",17),"Only the declared3019held-out anchors at192 enter population accuracy");
				final int n5=p.parseInt(10),n3=p.parseInt(11),r5=rawError(p.parseInt(6),p.parseInt(8),n5),r3=rawError(p.parseInt(7),p.parseInt(9),n3);
				require(p.termEquals(n5==0?"true":"false",12) && p.termEquals(n3==0?"true":"false",13)
					&& p.termEquals(n5==0 && n3==0?"true":"false",14),"Recorded exact flags must equal the native signed errors");
				final Assignment a=classify(m.taxid);scored++;m.assignment=a;m.raw5=r5;m.raw3=r3;m.net5=n5;m.net3=n3;m.scored=true;
				groups.get("GROUP\tALL\tALL\tALL").add(r5,r3,n5,n3);
				groups.get("GROUP\t"+a.group+"\tALL\tALL").add(r5,r3,n5,n3);
				final String key="PHYLUM\t"+a.group+'\t'+m.phylum+'\t'+a.phylum;Stats stats=groups.get(key);
				if(stats==null){stats=new Stats();groups.put(key,stats);}stats.add(r5,r3,n5,n3);
				detail.print(b.clear().append(m.id).tab().append(m.taxid).tab().append(m.phylum).tab().append(a.phylum).tab().append(a.group)
					.tab().append(a.ancestor==null?-1:a.ancestor.id).tab().append(a.ancestor==null?"NA":a.ancestor.name).tab().append(a.reason)
					.tab().append(r5).tab().append(r3).tab().append(n5).tab().append(n3).tab().append(a.lineage).nl());
			}
		}finally{final boolean error=in.close();RrnaResourceIO.finish(detail);require(!error,"Evaluation read failed");}
		require(seen==16110 && scored==3019,"Population report must preserve16110cohort rows and3019scored validation members: "+seen+"/"+scored);
		write(groups);writeOneBase(groups);if(gateFile!=null){writeGates();}System.err.println("MITO_ENDPOINT_POPULATION_REPORT_PASS cohort=16110 validation=3019 group_conservation=true scorer_not_caller=true");
	}
	private void loadRules(){
		final ByteFile in=RrnaResourceIO.open(ruleFile);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,"group\tancestor",ruleFile);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==2,"Every operational lineage rule needs one group and one exact ancestor name");final String group=p.parseString(0),name=p.parseString(1);
			require(group.equals("CHORDATE") || group.equals("ARTHROPOD") || group.equals("FUNGAL") || group.equals("ALGAL_LINEAGE"),"Unknown operational group");
			require(roots.put(name,new Root(group,name))==null,"Repeated configured ancestor: "+name);}
		}finally{require(!in.close(),"Lineage rules read failed");}
		require(roots.size()==18,"Expected three focal roots and15explicit algal roots");
		for(TaxNode n:tree.nodes){if(n==null){continue;}final Root r=roots.get(n.name);if(r!=null){require(r.node==null,"Ambiguous configured taxonomy name: "+r.name);r.node=n;}
			if(n.name.equals("Embryophyta")){require(landPlant==null,"Ambiguous land-plant exclusion");landPlant=n;}
			if(n.name.equals("Oomycota")){require(oomycete==null,"Ambiguous oomycete control");oomycete=n;}
			if(n.name.equals("Streptophyta")){require(broadGreen==null,"Ambiguous broad green-lineage control");broadGreen=n;}}
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".roots.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("group\tancestor_id\tancestor_name");for(String name:keys(roots)){final Root r=roots.get(name);require(r.node!=null,"Configured ancestor absent from frozen tree: "+name);out.print(b.clear().append(r.group).tab().append(r.node.id).tab().append(name).nl());}}
		finally{RrnaResourceIO.finish(out);}
	}
	private void loadEnrollment(){
		final ByteFile in=RrnaResourceIO.open(enrollment);final LineParser1 p=new LineParser1('\t');int panel=0;
		try{p.set(in.nextLine());require(p.terms()==19 && p.termEquals("role",0) && p.termEquals("phylum",4),"Expected accepted construction-ledger schema");
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==19,"Enrollment schema changed");if(p.termEquals("PANEL_INCLUDED",0)){panel++;continue;}
				require(p.termEquals("CONSTRUCTION",0),"Only construction members enter held-out population evaluation");final Member m=new Member(p.parseString(1),p.parseInt(3),p.parseString(4));
				require(m.taxid>0 && members.put(m.id,m)==null,"Duplicate/invalid original member identity");}
		}finally{require(!in.close(),"Enrollment read failed");}require(members.size()==16110 && panel==22,"Preserve construction16110 and separate panel22");
	}
	private Assignment classify(int taxid){
		// Unknown external hosts are an explicit report category. The nullable
		// lookup avoids TaxTree's newer-ID assertion; broken parent chains still fail.
		final Assignment a=new Assignment();TaxNode n=tree.getNode(taxid,true);final ByteBuilder lineage=new ByteBuilder();int depth=0;boolean land=false,broadGreen=false,ambiguous=false;
		if(n==null){a.reason="OTHER_UNRESOLVED_MISSING_HOST";a.lineage="NA";return a;}
		while(n!=null){require(++depth<1000,"Taxonomy parent cycle would invalidate lineage attribution for "+taxid);
			if(lineage.length()>0){lineage.append(';');}lineage.append(n.id).append(':').append(n.name);
			if(a.phylum.equals("RANK_ABSENT") && n.levelExtended==TaxTree.PHYLUM_E){a.phylum=n.id+":"+n.name;}
			land|=n.name.equals("Embryophyta");broadGreen|=n.name.equals("Streptophyta") || n.name.equals("Viridiplantae");
			final Root r=roots.get(n.name);if(r!=null && r.node.id==n.id){
				if(a.ancestor==null){a.group=r.group;a.ancestor=n;}else if(!a.group.equals(r.group)){ambiguous=true;}}
			if(n.id==n.pid){break;}n=tree.getNode(n.pid);require(n!=null,"Broken taxonomy parent chain for "+taxid);
		}
		require(!(land && a.group.equals("ALGAL_LINEAGE")),"The operational algal group must never include Embryophyta");
		if(ambiguous){a.group="OTHER";a.reason="OTHER_UNRESOLVED_MULTIPLE_GROUPS";}
		else if(a.ancestor!=null){a.reason="POSITIVE_CONFIGURED_ANCESTRY";}
		else if(land){a.reason="OTHER_LAND_PLANT";}else if(broadGreen){a.reason="OTHER_UNRESOLVED_GREEN_LINEAGE";}
		a.lineage=lineage.toString();return a;
	}
	/** Known taxonomy controls and independent coordinate examples run before any
	 * population numbers are emitted; failures cannot become plausible summaries. */
	private void policyFixture(){
		require(landPlant!=null && oomycete!=null && broadGreen!=null,"Frozen tree must resolve explicit negative lineage controls");
		require(classify(Integer.MAX_VALUE).reason.equals("OTHER_UNRESOLVED_MISSING_HOST"),"A newer/missing host must remain counted as unresolved");
		require(classify(landPlant.id).group.equals("OTHER") && classify(landPlant.id).reason.equals("OTHER_LAND_PLANT"),"Land plants must not be inferred to be algae");
		require(classify(broadGreen.id).reason.equals("OTHER_UNRESOLVED_GREEN_LINEAGE"),"Broad Streptophyta alone cannot prove algal ancestry");
		require(classify(roots.get("Charophyceae").node.id).group.equals("ALGAL_LINEAGE"),"A positively resolved charophyte branch must enter the operational algal bin");
		require(classify(oomycete.id).group.equals("OTHER") && classify(roots.get("Fungi").node.id).group.equals("FUNGAL"),"Oomycetes and true fungi must remain separate");
		require(rawError(40,70,20)==-10 && rawError(1040,1070,20)==-10,"Raw40/truth50/predicted70 means error-10 before and after a shared coordinate translation");
		final Stats fixture=new Stats();fixture.add(-1,2,0,-5);fixture.add(0,0,20,-20);
		require(java.util.Arrays.equals(fixture.within[0],new long[]{1,1,2,2,2,2}) && java.util.Arrays.equals(fixture.within[1],new long[]{0,0,0,1,1,2}),
			"Both-within must require both absolute errors, inclusive thresholds, and preserve the zero-tolerance exact case");
		System.err.println("MITO_POPULATION_POLICY_TEST_PASS missing_host_landplant_charophyte_oomycete_translation");
		System.err.println("MITO_WITHIN_TOLERANCE_TEST_PASS both_absolute_inclusive_zero_matches_exact");
	}
	private static int rawError(int raw,int predicted,int netError){return raw-predicted+netError;}
	/** Distinguish one wrong endpoint from two wrong endpoints within one base.
	 * The shift table is arithmetic on saved anchors, not a modified caller. */
	private void writeOneBase(ObjectMap<String,Stats> groups){
		final ByteStreamWriter pairs=RrnaResourceIO.writer(prefix+".onebase.tsv"),shift=RrnaResourceIO.writer(prefix+".shift3plus1.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			pairs.println("scope\tgroup\toriginal_phylum\tcurrent_phylum\tmembers\traw_error5\traw_error3\tpair_members");
			shift.println("scope\tgroup\toriginal_phylum\tcurrent_phylum\tmembers\traw_both_exact\tshifted_both_exact\texact_gained\texact_lost\traw_both_within10\tshifted_both_within10\twithin10_gained\twithin10_lost\tclaim");
			for(String key:keys(groups)){
				final Stats s=groups.get(key);long within=0;
				for(int a=-1;a<=1;a++){for(int z=-1;z<=1;z++){
					final long n=s.oneBase[(a+1)*3+z+1];within+=n;
					pairs.print(b.clear().append(key).tab().append(s.n).tab().append(a).tab().append(z).tab().append(n).nl());
				}}
				require(within==s.within[0][1] && s.oneBase[4]==s.exact[2],"Joint one-base bins must conserve the previously measured within-one and both-exact counts");
				pairs.print(b.clear().append(key).tab().append(s.n).append("\tOUTSIDE\tOUTSIDE\t").append(s.n-within).nl());
				shift.print(b.clear().append(key).tab().append(s.n).tab().append(s.exact[2]).tab().append(s.shiftExact).tab().append(s.shiftGainExact).tab().append(s.shiftLossExact)
					.tab().append(s.within[0][4]).tab().append(s.shiftNear).tab().append(s.shiftGainNear).tab().append(s.shiftLossNear).append("\tARITHMETIC_RAW_3PRIME_PLUS_ONE_NOT_CALLER").nl());
			}
		}finally{RrnaResourceIO.finish(pairs,shift);}
	}
	private void write(ObjectMap<String,Stats> groups){
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv"),offsets=RrnaResourceIO.writer(prefix+".offsets.tsv"),within=RrnaResourceIO.writer(prefix+".within.tsv");final ByteBuilder b=new ByteBuilder();long groupN=0,phylumN=0;
		try{summary.println("scope\tgroup\toriginal_phylum\tcurrent_phylum\tmembers\traw_exact5\traw_exact3\traw_both\tnet_exact5\tnet_exact3\tnet_both\tnet_both_fraction");
			offsets.println("scope\tgroup\toriginal_phylum\tcurrent_phylum\tstage\tend\tsigned_error\tmembers");
			within.println("scope\tgroup\toriginal_phylum\tcurrent_phylum\tmembers\ttolerance_bp\traw_both_within\tnet_both_within\traw5_within\traw3_within\tnet5_within\tnet3_within");
			for(String key:keys(groups)){final Stats s=groups.get(key);if(key.startsWith("PHYLUM\t")){phylumN+=s.n;}else if(!key.startsWith("GROUP\tALL\t")){groupN+=s.n;}
				b.clear().append(key).tab().append(s.n);for(long x:s.exact){b.tab().append(x);}b.tab();if(s.n==0){b.append("NA");}else{b.append(s.exact[5]/(double)s.n,8);}summary.print(b.nl());
				for(int i=0;i<TOLERANCES.length;i++){within.print(b.clear().append(key).tab().append(s.n).tab().append(TOLERANCES[i]).tab().append(s.within[0][i]).tab().append(s.within[1][i])
					.tab().append(s.endsWithin[0][0][i]).tab().append(s.endsWithin[0][1][i]).tab().append(s.endsWithin[1][0][i]).tab().append(s.endsWithin[1][1][i]).nl());}
				for(int j=0;j<4;j++){final IntList values=s.errors[j];require(values.size==s.n,"Every scored locus must contribute one signed error per end and stage");values.sort();
					for(int i=0;i<values.size;){final int value=values.get(i);int end=i+1;while(end<values.size && values.get(end)==value){end++;}
						offsets.print(b.clear().append(key).tab().append(j<2?"RAW":"NET").tab().append(j%2==0?"5prime":"3prime").tab().append(value).tab().append(end-i).nl());i=end;}}
			}
		}finally{RrnaResourceIO.finish(summary,offsets,within);}require(groupN==3019 && phylumN==3019,"Group and phylum partitions must independently conserve all3019validation members");
	}
	private void writeGates(){
		final ObjectMap<String,Gate> gates=new ObjectMap<String,Gate>(String.class,Gate.class);final ByteFile in=RrnaResourceIO.open(gateFile);final LineParser1 p=new LineParser1('\t');
		final ByteStreamWriter changes=RrnaResourceIO.writer(prefix+".gate_changes.tsv.gz");final ByteBuilder row=new ByteBuilder();
		changes.println("threshold\tlocus\thost\tgroup\traw_error5\traw_error3\tgated_error5\tgated_error3\tboth_exact_change\tboth_within10_change");
		try{RrnaResourceIO.header(in,"locus\thost\tthreshold\terror5\terror3\texact5\texact3\tboth_exact\tmoved5\tmoved3\tchosen_score5\tchosen_score3\tclaim",gateFile);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==13 && p.termEquals("SCORER_SEQUENTIAL_ABSOLUTE_SCORE_GATE_NOT_CALLER",12),"Gate rows must declare the sequential scorer policy");
				final Member m=members.get(p.parseString(0));require(m!=null && m.scored && m.taxid==p.parseInt(1),"Every gate row must join the same scored validation host/locus");
				final String threshold=p.parseString(2);final float value=Float.parseFloat(threshold);require(!Float.isNaN(value),"NaN cannot identify a gate");
				Gate g=gates.get(threshold);if(g==null){g=new Gate();gates.put(threshold,g);}require(g.seen.add(m.id),"Repeated locus at one threshold");
				final int d5=p.parseInt(3),d3=p.parseInt(4);
				require(p.termEquals(d5==0?"true":"false",5) && p.termEquals(d3==0?"true":"false",6) && p.termEquals(d5==0 && d3==0?"true":"false",7)
					&& p.termEquals(d5!=m.raw5?"true":"false",8) && p.termEquals(d3!=m.raw3?"true":"false",9)
					&& Float.isFinite(p.parseFloat(10)) && Float.isFinite(p.parseFloat(11)),"Gate exact/movement flags must reproduce signed errors and finite chosen scores");
				if(value==Float.NEGATIVE_INFINITY){require(d5==m.net5 && d3==m.net3,"Ungated cache must reproduce every original ungated endpoint");}
				if(value==Float.POSITIVE_INFINITY){require(d5==m.raw5 && d3==m.raw3,"All-reject gate must reproduce every raw anchor");}
				g.groups.get("ALL").add(m.raw5,m.raw3,d5,d3);g.groups.get(m.assignment.group).add(m.raw5,m.raw3,d5,d3);
				final String exact=change(m.raw5==0 && m.raw3==0,d5==0 && d3==0),near=change(near(m.raw5,m.raw3,10),near(d5,d3,10));
				if(!exact.equals("SAME") || !near.equals("SAME")){changes.print(row.clear().append(threshold).tab().append(m.id).tab().append(m.taxid).tab().append(m.assignment.group)
					.tab().append(m.raw5).tab().append(m.raw3).tab().append(d5).tab().append(d3).tab().append(exact).tab().append(near).nl());}
			}
		}finally{final boolean error=in.close();RrnaResourceIO.finish(changes);require(!error,"Gate rows read failed");}
		require(gates.get("-Infinity")!=null && gates.get("Infinity")!=null,"Both gate controls must be present in the saved experiment");
		final ArrayList<String> thresholds=keys(gates);Collections.sort(thresholds,(a,b)->Float.compare(Float.parseFloat(a),Float.parseFloat(b)));
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".gate_summary.tsv"),decisions=RrnaResourceIO.writer(prefix+".gate_decisions.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("threshold\tgroup\tmembers\traw_both_exact\tgated_both_exact\traw_both_within10\tgated_both_within10\tmoved_members\texact_nonregression\twithin10_nonregression\texact_gained\texact_lost\twithin10_gained\twithin10_lost");
			decisions.println("threshold\tall_five_lineages_nonregressed_exact_and_within10\tpositive_overall_gain\traw_both_exact\tgated_both_exact\traw_both_within10\tgated_both_within10\tmoved_members");
			for(String threshold:thresholds){final Gate g=gates.get(threshold);require(g.seen.size()==3019,"Every threshold must conserve all3019validation members");boolean safe=true;long sum=0;
				for(String group:keys(g.groups)){final Stats s=g.groups.get(group);final boolean exact=s.exact[5]>=s.exact[2],near=s.within[1][4]>=s.within[0][4];if(!group.equals("ALL")){safe&=exact && near;sum+=s.n;}
					out.print(b.clear().append(threshold).tab().append(group).tab().append(s.n).tab().append(s.exact[2]).tab().append(s.exact[5]).tab().append(s.within[0][4]).tab().append(s.within[1][4]).tab().append(s.moved).tab().append(exact).tab().append(near)
						.tab().append(s.gainExact).tab().append(s.lossExact).tab().append(s.gain10).tab().append(s.loss10).nl());}
				require(sum==3019,"Each threshold's five lineage bins must conserve the full population");final Stats s=g.groups.get("ALL");
				final boolean gain=s.moved>0 && (s.exact[5]>s.exact[2] || s.within[1][4]>s.within[0][4]);
				decisions.print(b.clear().append(threshold).tab().append(safe).tab().append(gain).tab().append(s.exact[2]).tab().append(s.exact[5]).tab().append(s.within[0][4]).tab().append(s.within[1][4]).tab().append(s.moved).nl());}
		}finally{RrnaResourceIO.finish(out,decisions);}
		System.err.println("MITO_ENDPOINT_GATE_REPORT_PASS thresholds="+gates.size()+" members_per_threshold=3019 raw_and_ungated_parity=true");
	}
	private static boolean near(int a,int b,int tolerance){return Math.abs((long)a)<=tolerance && Math.abs((long)b)<=tolerance;}
	private static String change(boolean before,boolean after){return before==after?"SAME":after?"GAIN":"LOSS";}
	private static <T> ArrayList<String> keys(ObjectMap<String,T> map){final ArrayList<String> list=new ArrayList<String>();for(String k:map.keys()){if(k!=null){list.add(k);}}Collections.sort(list);return list;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Root{Root(String g,String n){group=g;name=n;}final String group,name;TaxNode node;}
	private static final class Member{Member(String i,int t,String p){id=i;taxid=t;phylum=p;}final String id,phylum;final int taxid;boolean seen,scored;int raw5,raw3,net5,net3;Assignment assignment;}
	private static final class Assignment{String group="OTHER",phylum="RANK_ABSENT",reason="OTHER_NONFOCAL_LINEAGE",lineage;TaxNode ancestor;}
	private static final class Stats{
		void add(int r5,int r3,int n5,int n3){n++;if(r5==0){exact[0]++;}if(r3==0){exact[1]++;}if(r5==0 && r3==0){exact[2]++;}if(n5==0){exact[3]++;}if(n3==0){exact[4]++;}if(n5==0 && n3==0){exact[5]++;}errors[0].add(r5);errors[1].add(r3);errors[2].add(n5);errors[3].add(n3);
			final boolean rawExact=r5==0 && r3==0,netExact=n5==0 && n3==0,rawNear=near(r5,r3,10),netNear=near(n5,n3,10);
			if(near(r5,r3,1)){oneBase[(r5+1)*3+r3+1]++;}
			final boolean shiftedExact=r5==0 && r3==-1,shiftedNear=Math.abs((long)r5)<=10 && Math.abs((long)r3+1)<=10;
			if(shiftedExact){shiftExact++;}if(shiftedNear){shiftNear++;}
			if(!rawExact && shiftedExact){shiftGainExact++;}if(rawExact && !shiftedExact){shiftLossExact++;}
			if(!rawNear && shiftedNear){shiftGainNear++;}if(rawNear && !shiftedNear){shiftLossNear++;}
			if(!rawExact && netExact){gainExact++;}if(rawExact && !netExact){lossExact++;}if(!rawNear && netNear){gain10++;}if(rawNear && !netNear){loss10++;}
			if(r5!=n5 || r3!=n3){moved++;}for(int i=0;i<TOLERANCES.length;i++){final int t=TOLERANCES[i];if(near(r5,r3,t)){within[0][i]++;}if(near(n5,n3,t)){within[1][i]++;}
				if(Math.abs((long)r5)<=t){endsWithin[0][0][i]++;}if(Math.abs((long)r3)<=t){endsWithin[0][1][i]++;}if(Math.abs((long)n5)<=t){endsWithin[1][0][i]++;}if(Math.abs((long)n3)<=t){endsWithin[1][1][i]++;}}
		}
		long n,moved,gainExact,lossExact,gain10,loss10,shiftExact,shiftNear,shiftGainExact,shiftLossExact,shiftGainNear,shiftLossNear;
		final long[] oneBase=new long[9],exact=new long[6];final long[][] within=new long[2][TOLERANCES.length];final long[][][] endsWithin=new long[2][2][TOLERANCES.length];final IntList[] errors={new IntList(),new IntList(),new IntList(),new IntList()};
	}
	private static final class Gate{Gate(){for(String g:new String[]{"ALL","CHORDATE","ARTHROPOD","FUNGAL","ALGAL_LINEAGE","OTHER"}){groups.put(g,new Stats());}}final ObjectMap<String,Stats> groups=new ObjectMap<String,Stats>(String.class,Stats.class);final ObjectSet<String> seen=new ObjectSet<String>(String.class);}
	private String evaluation,enrollment,treeFile,ruleFile,prefix,gateFile;private TaxTree tree;private TaxNode landPlant,oomycete,broadGreen;
	private final ObjectMap<String,Root> roots=new ObjectMap<String,Root>(String.class,Root.class);
	private final ObjectMap<String,Member> members=new ObjectMap<String,Member>(String.class,Member.class);
	private static final String EVALUATION_HEADER="locus\thost\tsplit\tdisposition\tmodel_index\tradius\traw_start0\traw_stop0\tstart0\tstop0\terror5\terror3\texact5\texact3\tboth_exact\tstart_sites\tstop_sites\tclaim";
	private static final int[] TOLERANCES={0,1,2,5,10,20};
}
