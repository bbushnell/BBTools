package test.prok;

import prok.*;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Random;

/** Synthetic plumbing and both-strand anchor fixtures, not learned 18S resources.
 * @author Raiden
 */
public final class Euk18sRuntimeConfigTest {
	public static void main(String[] args)throws Exception{
		final Path dir=Files.createTempDirectory("euk18s_config_");
		final byte[] gene=new byte[1803],alphabet=bytes("ACGT");final Random random=new Random(90317);
		for(int i=0;i<gene.length;i++){gene[i]=alphabet[random.nextInt(4)];}
		write(dir,"refs.fa",">pilot\n"+new String(gene,StandardCharsets.US_ASCII)+"\n");
		write(dir,"seeds.fa",">one\n"+new String(gene,200,17,StandardCharsets.US_ASCII)+"\n>two\n"+new String(gene,1300,17,StandardCharsets.US_ASCII)+"\n");
		write(dir,"cutoffs.tsv","family\tmodel\tidpass\tidborderline\neuk18S\tpilot\t0.8\t0.8\n");
		write(dir,"vote_gene.fa",">synthetic\nAAAAA"+new String(gene,StandardCharsets.US_ASCII)+"CCCCCCC\n");
		SeedOffsetTableBuilder.build(dir.resolve("vote_gene.fa").toString(),dir.resolve("seeds.fa").toString(),dir.resolve("offsets.tsv").toString(),"00000000000000000000",17,true);
		final ArrayList<NcrnaFamily> empty=new ArrayList<NcrnaFamily>();new Euk18sRuntimeConfig().apply(empty);
		check(empty.isEmpty(),"Disabled family must perform no loading or registration");
		final Euk18sRuntimeConfig off=new Euk18sRuntimeConfig();off.parse("euk18sconsensus",dir.resolve("missing").toString());
		reject(()->off.apply(empty),"require euk18s=t");
		final NcrnaFamily alignment=load(dir,false),voted=load(dir,true);
		check(alignment.kLong==17&&alignment.indexK==9&&alignment.seedMinHits==2&&!alignment.seedDistinct&&alignment.windowPad==1903
			&&alignment.outputType==ProkObject.r18S&&alignment.modelThresholds.pass(0)==.8f,"Registration must bind actual resource controls");
		final byte[] contig=new byte[2203];for(int i=0;i<contig.length;i++){contig[i]=alphabet[random.nextInt(4)];}
		System.arraycopy(gene,0,contig,200,gene.length);
		for(int strand=0;strand<2;strand++){
			final NcrnaScavenger a=RrnaEndpointCallerFeaturesTest.scavenger(alignment),v=RrnaEndpointCallerFeaturesTest.scavenger(voted);
			final ArrayList<Orf> ca=a.scavenge("fixture",contig,strand,new ArrayList<int[]>());
			final ArrayList<Orf> cv=v.scavenge("fixture",contig,strand,new ArrayList<int[]>());
			check(ca.size()==1&&cv.size()==1,"Both raw endpoint policies must detect the exact synthetic model");
			check(ca.get(0).start==200&&ca.get(0).stop==200+gene.length-1,"Alignment mode must retain Quantum endpoints");
			check(cv.get(0).start==195&&cv.get(0).stop==200+gene.length+6,"Voting mode must use bound offsets in strand-oriented coordinates");
			check(v.votedStarts==1&&v.votedStops==1,"The fixture must actually execute both endpoint votes");
		}
		// Observation must see voted anchors, not the preceding alignment anchors.
		final RrnaPositionalKmerTable.Table five=new RrnaPositionalKmerTable.Table("pilot","5prime",6,0,.5,16);
		final RrnaPositionalKmerTable.Table three=new RrnaPositionalKmerTable.Table("pilot","3prime",9,-8,.5,16);
		five.add(contig,200);three.add(contig,200+gene.length-1);five.finish();three.finish();final int[] observed={0};
		voted.setRrnaEndpointFeatures(new RrnaPositionalKmerTable.Table[]{five},new RrnaPositionalKmerTable.Table[]{three},
			(c,s,m,start,stop,f,fu,t,tu)->{check(start==195&&stop==200+gene.length+6&&f.length==36&&t.length==36,"Resource-sized features must observe the selected raw votes");observed[0]++;});
		RrnaEndpointCallerFeaturesTest.scavenger(voted).scavenge("fixture",contig,0,new ArrayList<int[]>());
		check(observed[0]==1,"Expected one actual voted feature capture");
		endpointDirectoryFixture(dir,contig);
		final NcrnaFamily malformed=load(dir,false);malformed.voteWindows=true;
		malformed.setRrnaEndpointFeatures(new RrnaPositionalKmerTable.Table[]{five},new RrnaPositionalKmerTable.Table[]{three},voted.rrnaEndpointFeatureSink);
		reject(()->RrnaEndpointCallerFeaturesTest.scavenger(malformed),"optional bound votes");
		final Euk18sRuntimeConfig duplicate=config(dir,false);final ArrayList<NcrnaFamily> registered=new ArrayList<NcrnaFamily>();registered.add(alignment);
		reject(()->duplicate.apply(registered),"Duplicate 18S");check(registered.size()==1,"Duplicate failure must preserve the original family");
		final Euk18sRuntimeConfig missing=config(dir,false);missing.parse("euk18smodelcutoffs",dir.resolve("absent.tsv").toString());
		reject(()->missing.apply(empty),"Missing required euk18S");check(empty.isEmpty(),"Resource failure must not publish a partial family");
		final Euk18sRuntimeConfig noVotes=config(dir,false);noVotes.parse("euk18srawends","vote");noVotes.parse("euk18svoteoffsets",dir.resolve("absent.offsets").toString());reject(()->noVotes.apply(empty),"Missing required euk18S");
		final Euk18sRuntimeConfig wrongSeed=config(dir,true);wrongSeed.parse("euk18svoteseedsha80","11111111111111111111");reject(()->wrongSeed.apply(empty),"sha80 mismatch");
		reject(()->config(dir,false).parse("euk18sseedminhits","0"),"must be >=1");
		final Euk18sRuntimeConfig seedMode=config(dir,false);
		seedMode.parse("euk18sseedminhits","1");seedMode.parse("euk18sseeddistinct","t");
		final ArrayList<NcrnaFamily> seedFamilies=new ArrayList<NcrnaFamily>();seedMode.apply(seedFamilies);
		final NcrnaScavenger seedCaller=RrnaEndpointCallerFeaturesTest.scavenger(seedFamilies.get(0));
		check(seedFamilies.get(0).seedMinHits==1 && seedFamilies.get(0).seedDistinct
			&& seedCaller.minKmerHits==1 && seedCaller.seedDistinct,"Minimum1 and distinct mode must reach the actual caller");
		reject(()->config(dir,false).parse("euk18sindexmargin","-2"),"must be >=-1");
		final Euk18sRuntimeConfig conflict=config(dir,false);conflict.parse("euk18stopn","3");conflict.parse("euk18sindexmargin","0");
		reject(()->conflict.apply(empty),"Choose euk18sindexmargin or euk18stopn");
		marginFixture(dir,gene,contig);
		fractionFixture(dir,contig);
		voteWindowFixture(dir,contig);
		// The independent registration enters exactly the same shared K17 index.
		GeneCaller.ncrnaFamilies.clear();GeneCaller.ncrnaFamilies.add(alignment);ProkObject.calltRNA=false;
		GeneCaller.initializeConservedRnaSeedIndex();check(GeneCaller.conservedRnaSeedSlotCount()==1,"18S must share the generic K17 scanner");
		System.out.println("PASS Euk18sRuntimeConfigTest registration shared_scan both_raw_modes both_strands voted_features resource_binding disabled_path");
	}
	/** Exercise the optional CLI directories against the preceding direct-loader API.
	 * Random nets are plumbing fixtures, not calibrated or packaged18S endpoints. */
	private static void endpointDirectoryFixture(Path dir,byte[] contig)throws Exception{
		assert(contig!=null && contig.length>2003) : "Synthetic endpoint tables need both full gene boundaries";
		if(ml.Function.TYPE_RATES_CUM==null){ml.Function.normalizeTypeRates();}
		final Path root=Files.createDirectory(dir.resolve("endpoint_directories"));
		Files.createDirectory(root.resolve("networks"));
		final structures.ByteBuilder tables=new structures.ByteBuilder().append(RrnaEndpointResourceLoader.HEADER).nl();
		final structures.ByteBuilder nets=new structures.ByteBuilder().append("model\tend\tseed\tfile\tsha80").nl();
		for(int end=0; end<2; end++){
			final String label=end==0 ? "5prime" : "3prime";
			final int k=end==0 ? 6 : 9, offset=end==0 ? -3 : -8;
			final RrnaPositionalKmerTable.Table table=new RrnaPositionalKmerTable.Table("pilot",label,k,offset,.5,16);
			table.add(contig,end==0 ? 200 : 2002);table.finish();
			table.write(root.resolve("pilot."+label+".tsv.gz").toString());
			// Existing synthetic loader fixtures use an explicit all-zero provenance placeholder.
			tables.append("pilot").tab().append(label).tab().append(k).tab().append(16).tab().append(offset).append("\t0.5\t1\t00000000000000000000").nl();
			final ml.CellNet net=new ml.CellNet(new int[]{36,33},74003,1f,1f,1,new ArrayList<String>());net.randomize();
			final String file="networks/pilot."+label+".s74003.bbnet";
			final fileIO.ByteStreamWriter writer=RrnaResourceIO.writer(root.resolve(file).toString());
			writer.print(net.toBytes());RrnaResourceIO.finish(writer);
			nets.append("pilot").tab().append(label).append("\t74003\t").append(file).append("\t00000000000000000000").nl();
		}
		final fileIO.ByteStreamWriter tw=RrnaResourceIO.writer(root.resolve("TABLES.tsv").toString());
		tw.print(tables);RrnaResourceIO.finish(tw);
		final fileIO.ByteStreamWriter nw=RrnaResourceIO.writer(root.resolve("MANIFEST.tsv").toString());
		nw.print(nets);RrnaResourceIO.finish(nw);
		final Euk18sRuntimeConfig c=config(dir,false);c.parse("euk18sendpoint","t");
		c.parse("euk18sendpointtables",root.toString());c.parse("euk18sendpointnets",root.toString());
		final ArrayList<NcrnaFamily> list=new ArrayList<NcrnaFamily>();c.apply(list);
		check(list.size()==1,"Explicit endpoint resources must register one complete18S family");
		final NcrnaFamily family=list.get(0);
		final RrnaEndpointCallerFeatures.Resources actual=family.rrnaEndpointFeatures;
		final RrnaEndpointCallerFeatures.Resources old=RrnaEndpointNetworkLoader.load(root.toString(),
			RrnaEndpointResourceLoader.load("euk18S",root.toString(),family.modelNames,family.library));
		check(actual!=null && actual.radius==16 && actual.inputs==36 && !actual.acceptedVerifierIdentity,
			"External18S endpoints retain candidate-span identity and their resource-defined geometry");
		final float[][] observed=new float[2][];
		family.rrnaEndpointFeatureSink=(name,strand,m,start,stop,f,fu,t,tu)->{observed[0]=f.clone();observed[1]=t.clone();};
		// The worker copies family endpoint state at construction.
		final NcrnaScavenger reporting=RrnaEndpointCallerFeaturesTest.scavenger(family);
		final Orf candidate=new Orf("external_endpoint",203,1998,0,0,contig,false,ProkObject.r18S);
		candidate.ncrnaFamily="euk18S";candidate.trnaModel="pilot";
		reporting.captureRrnaEndpointFeatures(candidate,contig,0);
		check(observed[0]!=null && observed[1]!=null
			&& candidate.start==203+RrnaEndpointCallerFeatures.predict(old.fiveNets[0],observed[0])-16
			&& candidate.stop==1998+RrnaEndpointCallerFeatures.predict(old.threeNets[0],observed[1])-16,
			"Actual migrated runtime inference must retain the prior directory-loader endpoint predictions");
		final Euk18sRuntimeConfig missing=config(dir,false);missing.parse("euk18sendpoint","t");
		reject(()->missing.apply(new ArrayList<NcrnaFamily>()),"explicit table and network directories");
	}
	private static void voteWindowFixture(Path dir,byte[] contig){
		final Euk18sRuntimeConfig c=config(dir,true);c.parse("euk18srawends","alignment");
		c.parse("euk18svotewindows","t");c.parse("euk18svotebeforepadding","t");c.parse("euk18swindowpad","3109");
		final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();c.apply(families);final NcrnaFamily f=families.get(0);
		check(f.windowPad==3109 && f.voteBeforePadding && f.voteWindows && !f.voteEnds,"Explicit padding/window mode must reach the family without changing raw endpoints");
		for(int strand=0;strand<2;strand++){
			final NcrnaScavenger s=RrnaEndpointCallerFeaturesTest.scavenger(f);
			check(s.voteBeforePadding,"GeneCaller control bridge must propagate the new window policy");
			final ArrayList<Orf> calls=s.scavenge("vote_first",contig,strand,new ArrayList<int[]>());
			check(calls.size()==1 && calls.get(0).start==200 && calls.get(0).stop==200+1803-1 && s.votedWindows==1,
				"Actual opted-in caller must execute voted windows while retaining alignment endpoints on each oriented strand");
		}
		final Euk18sRuntimeConfig missing=config(dir,true);missing.parse("euk18svotebeforepadding","t");
		reject(()->missing.apply(new ArrayList<NcrnaFamily>()),"explicit measured windowpad");
		reject(()->config(dir,false).parse("euk18swindowpad","-1"),"must be >=0");
		check(!load(dir,false).voteBeforePadding,"Omitted new flag must retain the legacy window builder");
	}
	/** The repetitive decoy ranks first by occurrence-weighted k9 score but fails
	 * alignment. A wide gap admits the real model second and exits before its twin. */
	private static void marginFixture(Path dir,byte[] gene,byte[] contig)throws Exception{
		final String sequence=new String(gene,StandardCharsets.US_ASCII);
		write(dir,"margin_refs.fa",">decoy\n"+sequence+sequence+"\n>pilot\n"+sequence+"\n>twin\n"+sequence+"\n");
		write(dir,"margin_cutoffs.tsv","family\tmodel\tidpass\tidborderline\neuk18S\tdecoy\t0.99\t0.99\neuk18S\tpilot\t0.99\t0.99\neuk18S\ttwin\t0.99\t0.99\n");
		for(int margin:new int[]{-1,0,Integer.MAX_VALUE}){
			final Euk18sRuntimeConfig c=config(dir,false);
			c.parse("euk18sconsensus",dir.resolve("margin_refs.fa").toString());
			c.parse("euk18smodelcutoffs",dir.resolve("margin_cutoffs.tsv").toString());
			c.parse("euk18sindexmargin",Integer.toString(margin));
			final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();c.apply(families);
			final NcrnaScavenger s=RrnaEndpointCallerFeaturesTest.scavenger(families.get(0));
			check(s.indexScoreMargin==margin,"CLI margin must reach the production scavenger");
			final int[] attempts={0};
			s.setModelAttemptSink((name,strand,pass,start,stop,rank,model,shared,id,length,accepted)->{
				check(model==attempts[0] && rank==model+1,"Decoy must rank first, then pilot; twin must never align");attempts[0]++;
			});
			final ArrayList<Orf> calls=s.scavenge("margin",contig,0,new ArrayList<int[]>());
			check(calls.size()==(margin==0?0:1) && attempts[0]==(margin==0?1:2),
				"Margin must change actual alignments and preserve first-accepted early exit; margin="+margin+" calls="+calls.size()+" attempts="+attempts[0]);
		}
	}
	/** Actual config->family->scavenger->index execution, using the margin fixture's
	 * doubled decoy and two identical true models. The high B excludes the latter. */
	private static void fractionFixture(Path dir,byte[] contig){
		for(int margin:new int[]{-1,Integer.MAX_VALUE}){
			String omitted=null;
			for(String fraction:new String[]{null,"0","0.9"}){
				final Euk18sRuntimeConfig c=fractionConfig(dir,margin,fraction,96);final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();c.apply(families);
				final NcrnaFamily f=families.get(0);final NcrnaScavenger s=RrnaEndpointCallerFeaturesTest.scavenger(f);final boolean high="0.9".equals(fraction);
				check(f.indexTopN==f.library.length && s.indexScoreMargin==margin && s.adaptiveMinHits==high
					&& s.adaptQFrac==0 && s.indexMinHitsDefault==96 && (!high || (s.adaptFloor==96 && s.adaptTopFrac==.9f)),
					"A/B/F controls must reach the actual caller without query-fraction filtering or a hidden model cap");
				final StringBuilder events=new StringBuilder();final int[] attempts={0};
				s.setModelAttemptSink((name,strand,pass,start,stop,rank,model,shared,id,length,accepted)->{
					events.append(rank).append(':').append(model).append(':').append(shared).append(':').append(id).append(':').append(length).append(':').append(accepted).append(';');attempts[0]++;
				});
				final ArrayList<Orf> calls=s.scavenge("fraction",contig,0,new ArrayList<int[]>());
				check(calls.size()==(high?0:1) && attempts[0]==(high?1:2),"Nonzero B must prune actual downstream model attempts with A both off and enabled");
				for(Orf call:calls){events.append(call.start).append(':').append(call.stop).append(':').append(call.trnaModel).append(':').append(call.orfScore);}
				if(fraction==null){omitted=events.toString();}else if(!high){check(omitted.equals(events.toString()),"Omitted B and explicit zero must preserve exact model attempts, identities and called endpoints/scores");}
			}
			for(String fraction:new String[]{"0","0.9"}){
				final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();fractionConfig(dir,margin,fraction,Integer.MAX_VALUE).apply(families);
				final NcrnaScavenger s=RrnaEndpointCallerFeaturesTest.scavenger(families.get(0));
				check(s.scavenge("floor",contig,0,new ArrayList<int[]>()).isEmpty() && s.alignmentCount()==0,
					"Absolute F must reject below-floor windows before alignment even with A off and keep-best index fallback");
			}
		}
		for(String invalid:new String[]{"-0.1","1.01","NaN","Infinity"}){reject(()->config(dir,false).parse("euk18sindexfrac",invalid),"outside [0.0,1.0]");}
		final Euk18sRuntimeConfig capped=fractionConfig(dir,-1,"0.9",96);capped.parse("euk18stopn","1");reject(()->capped.apply(new ArrayList<NcrnaFamily>()),"do not specify euk18stopn");
		final Euk18sRuntimeConfig exact=config(dir,false);exact.parse("euk18sindexfrac","1");
	}
	private static Euk18sRuntimeConfig fractionConfig(Path dir,int margin,String fraction,int floor){
		final Euk18sRuntimeConfig c=config(dir,false);c.parse("euk18sconsensus",dir.resolve("margin_refs.fa").toString());c.parse("euk18smodelcutoffs",dir.resolve("margin_cutoffs.tsv").toString());
		c.parse("euk18sindexmargin",Integer.toString(margin));c.parse("euk18sindexminhits",Integer.toString(floor));if(fraction!=null){c.parse("euk18sindexfrac",fraction);}return c;
	}
	static Euk18sRuntimeConfig config(Path dir,boolean vote){
		final Euk18sRuntimeConfig c=new Euk18sRuntimeConfig();c.parse("EuK18S","t");
		c.parse("euk18sjoined","f");
		// These historical fixtures isolate explicit Quantum/window policies. Plain
		// euk18s=t is separately checked against the frozen accepted shred recipe.
		c.parse("euk18saligner","quantum");c.parse("euk18smodelclip","f");
		c.parse("euk18svotewindows","f");c.parse("euk18svotebeforepadding","f");
		c.parse("euk18sseeddistinct","f");c.parse("euk18swindowpad","1903");
		c.parse("euk18sindexminhits","1");c.parse("euk18sindexmargin","-1");c.parse("euk18sindexfrac","0");
		c.parse("euk18sconsensus",dir.resolve("refs.fa").toString());c.parse("euk18sseeds",dir.resolve("seeds.fa").toString());c.parse("euk18smodelcutoffs",dir.resolve("cutoffs.tsv").toString());
		if(vote){c.parse("euk18srawends","vote");c.parse("euk18svoteoffsets",dir.resolve("offsets.tsv").toString());c.parse("euk18svoteseedsha80","00000000000000000000");}return c;
	}
	static NcrnaFamily load(Path dir,boolean vote){final ArrayList<NcrnaFamily> f=new ArrayList<NcrnaFamily>();config(dir,vote).apply(f);check(f.size()==1,"One complete family must register");return f.get(0);}
	static void write(Path root,String name,String data)throws Exception{Files.write(root.resolve(name),bytes(data));}
	static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	static void reject(Runnable run,String why){try{run.run();}catch(IllegalArgumentException e){check(e.getMessage().contains(why),"Unexpected failure: "+e);return;}throw new AssertionError("Expected "+why);}
}
