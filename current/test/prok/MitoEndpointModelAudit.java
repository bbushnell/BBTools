package test.prok;

import prok.MitoSsuModelPolicy;
import prok.RrnaResourceIO;
import prok.TrnaConsensusBuilder;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntHashMap;
import template.ThreadWaiter;

/** Saved full-cohort model offsets and bounded, fixed-assignment terminal edits.
 * ByteFile-template parsing, native FASTA streaming and independent worker batches.
 * No caller, model reassignment, resource mutation, or learned boundary network.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden */
public final class MitoEndpointModelAudit {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, MitoEndpointModelAudit.class, false);
		try{final MitoEndpointModelAudit x=new MitoEndpointModelAudit();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		final Parser p=new Parser();
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg, '=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("mode")){mode=v;}else if(k.equals("enrollment")){enrollment=v;}
			else if(k.equals("consensus")){consensus=v;}else if(k.equals("corpus")){corpus=v;}
			else if(k.equals("edits")){edits=v;}else if(k.equals("outprefix")){prefix=v;}
			else if(k.equals("model")){supportModel=Integer.parseInt(v);}else if(k.equals("end")){supportEnd=Integer.parseInt(v);}
			else if(k.equals("expected")){expectedSupport=Integer.parseInt(v);}
			else if(!p.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		assignments=p.in1;require(assignments!=null && enrollment!=null && consensus!=null && prefix!=null && !p.append && p.maxReads<0,
			"Use complete frozen assignments, enrollment, unchanged library and fresh outputs");
		require(mode.equals("counts") || mode.equals("support") || mode.equals("test"),"Declare counts, unedited source support, or bounded fixed-assignment test");
		require(mode.equals("counts")?edits==null && corpus==null:mode.equals("support")?edits==null && corpus!=null:edits!=null && corpus!=null,"Support mode never accepts an edit roster");
		require(!mode.equals("support") || supportModel>=5 && supportModel<158 && (supportEnd==3 || supportEnd==5) && expectedSupport>0,"Source support needs an explicit mitochondrial model/end and census denominator");
		Shared.TRIM_READ_DESCRIPTION=false;stream.Read.TO_UPPER_CASE=true;
		ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;ReadWrite.ZIPLEVEL=9;ReadWrite.setZipThreads(1);
	}
	private void run(){
		library=TrnaConsensusBuilder.loadLibrary(consensus);names=TrnaConsensusBuilder.lastLoadedNames;
		require(library.length==158 && MitoSsuModelPolicy.coverage("16S",names),"Keep the original five plus all153shipped mitochondrial models");
		readEnrollment();readAssignments();
		if(mode.equals("counts")){writeCounts();return;}
		if(mode.equals("support")){
			RrnaResourceIO.fresh(prefix,".support_members.tsv.gz",".support_summary.tsv");
			variants=new byte[158][][];variants[supportModel]=new byte[][]{library[supportModel].clone()};variantNames[supportModel]=new String[]{"BASELINE"};ends[supportModel]=supportEnd;
		}else{RrnaResourceIO.fresh(prefix,".test_members.tsv.gz",".test_summary.tsv",".test_bases.tsv",".control.tsv");readEdits();}
		loadWindows();
		final AtomicInteger next=new AtomicInteger();final ArrayList<Worker> workers=new ArrayList<Worker>();
		for(int i=0;i<Shared.threads();i++){workers.add(new Worker(next));}
		require(ThreadWaiter.startAndWait(workers),"All read-only alignment workers must terminate before reporting");
		for(Worker w:workers){if(w.failure!=null){throw new IllegalStateException("Model diagnostic failed",w.failure);}require(w.success,"Incomplete worker");}
		if(mode.equals("support")){writeSupport();}else{writeTests();}
	}
	private void readEnrollment(){
		final ByteFile in=RrnaResourceIO.open(enrollment);final LineParser1 p=new LineParser1('\t');int panel=0;
		try{p.set(in.nextLine());require(p.terms()==19 && p.termEquals("role",0) && p.termEquals("phylum",4),"Expected accepted coverage enrollment schema");
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==19,"Enrollment width changed");if(p.termEquals("PANEL_INCLUDED",0)){panel++;continue;}
				require(p.termEquals("CONSTRUCTION",0),"Never silently include evaluation-panel members in construction counts");
				final Row r=new Row();r.id=p.parseString(1);r.host=p.parseInt(3);r.phylum=p.parseString(4);
				require(r.host>0 && byId.put(r.id,r)==null,"Construction membership must retain unique IDs and positive hosts");
			}
		}finally{require(!in.close(),"Enrollment read failed");}
		require(byId.size()==16110 && panel==22,"Preserve16110construction members and separate22panel rows");
	}
	private void readAssignments(){
		final ByteFile in=RrnaResourceIO.open(assignments);final LineParser1 p=new LineParser1('\t');int usable=0,validation=0;
		try{RrnaResourceIO.header(in,HEADER,assignments);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==21,"Saved assignment schema changed");final Row r=byId.get(p.parseString(0));
				require(r!=null && !r.seen && p.parseInt(2)==r.host,"Assignments must conserve every enrolled source once");r.seen=true;rows.add(r);
				require(p.termEquals("TRAIN",3) || p.termEquals("VALIDATION",3),"Every source retains its declared host split");r.split=p.termEquals("TRAIN",3)?1:2;
				r.disposition=p.parseString(6);r.model=p.parseInt(7);r.ts=p.parseInt(13);r.te=p.parseInt(14);r.s=p.parseInt(15);r.e=p.parseInt(16);
				r.eligible=r.disposition.equals("COMPLETE_SOURCE_ASSOCIATED");
				if(r.eligible){
					require(r.model>=5 && r.model<library.length && p.termEquals(names[r.model],8) && p.termEquals("false",19) && p.termEquals("false",20)
						&& r.s>=0 && r.e>=r.s && r.ts>=0 && r.te>=r.ts,"Eligible anchors must retain full source bounds and exact shipped model identity");
					require(p.parseInt(17)==r.ts-r.s && p.parseInt(18)==r.te-r.e,"Saved corrections must equal truth minus raw, not signed errors");
					r.identity=p.parseFloat(12);require(Float.isFinite(r.identity),"Saved genomic identity must be finite");usable++;if(r.split==2){validation++;}
				}
			}
		}finally{require(!in.close(),"Assignment read failed");}
		require(rows.size()==16110 && usable==15084 && validation==3019,"Complete population15084=12065TRAIN+3019VALIDATION must be preserved");
	}
	private void writeCounts(){
		RrnaResourceIO.fresh(prefix,".models.tsv",".phyla.tsv",".dispositions.tsv",".offsets.tsv",".modes.tsv");
		final Counts[][] all=new Counts[158][3];for(int m=5;m<158;m++){for(int split=0;split<3;split++){all[m][split]=new Counts();}}
		final ObjectMap<String,Counts> phyla=new ObjectMap<String,Counts>(String.class,Counts.class),states=new ObjectMap<String,Counts>(String.class,Counts.class);
		for(Row r:rows){for(int split:new int[]{0,r.split}){final String key=SPLITS[split]+'\t'+r.disposition;Counts c=states.get(key);if(c==null){c=new Counts();states.put(key,c);}c.n++;}
			if(!r.eligible){continue;}final int a=r.s-r.ts,z=r.e-r.te;
			for(int split:new int[]{0,r.split}){all[r.model][split].add(a,z);final String key=r.model+"\t"+names[r.model]+'\t'+SPLITS[split]+'\t'+r.phylum;
				Counts c=phyla.get(key);if(c==null){c=new Counts();phyla.put(key,c);}c.add(a,z);}
		}
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".models.tsv"),ph=RrnaResourceIO.writer(prefix+".phyla.tsv"),ds=RrnaResourceIO.writer(prefix+".dispositions.tsv"),hist=RrnaResourceIO.writer(prefix+".offsets.tsv"),modes=RrnaResourceIO.writer(prefix+".modes.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("model_index\tmodel\tsplit\tmodel_length\t"+COUNTS);ph.println("model_index\tmodel\tsplit\tinput_phylum\t"+COUNTS);ds.println("split\tdisposition\tmembers");
			hist.println("model_index\tmodel\tsplit\tend\tsigned_error\tmembers");
			modes.println("model_index\tmodel\tsplit\tend\tmembers\tmodal_error\tmodal_count\tmodal_ties\tmodal_fraction");
			for(int m=5;m<158;m++){for(int split=0;split<3;split++){b.clear().append(m).tab().append(names[m]).tab().append(SPLITS[split]).tab().append(library[m].length);all[m][split].append(b);out.print(b.nl());}}
			for(int m=5;m<158;m++){for(int split=0;split<3;split++){final Counts c=all[m][split];for(int end=0;end<2;end++){
				final IntHashMap h=c.errors[end];final int[] values=h.toArray();Arrays.sort(values);int max=0,best=0,ties=0;long total=0;
				for(int value:values){final int n=h.get(value);total+=n;if(n>max){max=n;best=value;ties=1;}else if(n==max){ties++;}
					hist.print(b.clear().append(m).tab().append(names[m]).tab().append(SPLITS[split]).tab().append(end==0?5:3).tab().append(value).tab().append(n).nl());}
				require(total==c.n,"Every endpoint histogram must conserve the complete per-model/split membership");
				b.clear().append(m).tab().append(names[m]).tab().append(SPLITS[split]).tab().append(end==0?5:3).tab().append(c.n).tab();
				if(c.n==0){b.append("NA\t0\t0\tNA");}else{b.append(best).tab().append(max).tab().append(ties).tab().append(max/(double)c.n,9);}modes.print(b.nl());
			}}}
			for(String key:keys(phyla)){b.clear().append(key);phyla.get(key).append(b);ph.print(b.nl());}
			for(String key:keys(states)){ds.print(b.clear().append(key).tab().append(states.get(key).n).nl());}
		}finally{RrnaResourceIO.finish(out,ph,ds,hist,modes);}
		System.err.println("MITO_MODEL_OFFSET_COUNTS_PASS cohort=16110 eligible=15084 train=12065 validation=3019 models=153 alignments=0");
	}
	private void readEdits(){
		variants=new byte[158][][];final ByteFile in=RrnaResourceIO.open(edits);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,"model_index\tmodel\tend\ttrim\textend_bases",edits);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==5,"Edits require fixed identity/end and declared source-supported alternatives");final int m=p.parseInt(0),end=p.parseInt(2),trim=p.parseInt(3);
				require(m>=5 && m<158 && p.termEquals(names[m],1) && variants[m]==null && (end==3 || end==5),"No original-model edit, duplicate model or implicit endpoint");
				final String bases=p.termEquals(".",4)?"":p.parseString(4);require(trim>=0 && trim<=2 && (trim>0 || !bases.isEmpty()),"Only declared one/two-base trims or explicitly supported single-base extensions");
				final int count=1+(trim>0?1:0)+bases.length();variants[m]=new byte[count][];variantNames[m]=new String[count];ends[m]=end;
				variants[m][0]=library[m].clone();variantNames[m][0]="BASELINE";int at=1,mask=0;
				if(trim>0){variants[m][at]=terminalEdit(library[m],end,trim,(byte)0);variantNames[m][at++]="TRIM_"+trim;}
				for(int i=0;i<bases.length();i++){final byte base=(byte)bases.charAt(i);final int bit="ACGTN".indexOf(base);require(bit>=0 && (mask&(1<<bit))==0,"Extensions must be unique explicit A/C/G/T or owner-requested neutral N");mask|=1<<bit;
					if(base!='N'){literalComparator[m]=at;}
					variants[m][at]=terminalEdit(library[m],end,0,base);variantNames[m][at++]="EXTEND_"+(char)base;}
				if(bases.indexOf('N')>=0){require(trim==0 && bases.length()==2 && literalComparator[m]>0,"Neutral-N experiment requires exactly one declared literal comparator plus N, with no trim");}
				selected++;
			}
		}finally{require(!in.close(),"Edit roster read failed");}require(selected>=1 && selected<=3,"Bounded diagnostic permits only one to three explicitly chosen models");
	}
	/** Baseline, trim one, extend A/C/G/T. No length or base is silently optimized. */
	public static byte[][] terminalVariants(byte[] model,int end){
		require(model!=null && model.length>16 && (end==3 || end==5),"Terminal edits need an intact model and explicit5'or3'end");
		final byte[][] v=new byte[6][];v[0]=model.clone();v[1]=Arrays.copyOfRange(model,end==5?1:0,model.length-(end==3?1:0));
		for(int i=0;i<4;i++){v[i+2]=new byte[model.length+1];System.arraycopy(model,0,v[i+2],end==5?1:0,model.length);v[i+2][end==5?0:model.length]=BASES[i];}
		return v;
	}
	/** Copies one declared terminal edit; never edits the packaged library in place. */
	public static byte[] terminalEdit(byte[] model,int end,int trim,byte extend){
		require(model!=null && (end==3 || end==5) && trim>=0 && trim<=2 && model.length-trim>16
			&& (trim>0?extend==0:"ACGTN".indexOf(extend)>=0),"Choose one small trim OR an explicitly declared canonical/neutral-N extension");
		if(trim>0){return Arrays.copyOfRange(model,end==5?trim:0,model.length-(end==3?trim:0));}
		final byte[] result=new byte[model.length+1];System.arraycopy(model,0,result,end==5?1:0,model.length);result[end==5?0:model.length]=extend;return result;
	}
	private void loadWindows(){
		RrnaResourceIO.stream(corpus,read->{final Row r=byId.get(RrnaResourceIO.first(read.id));require(r!=null && !r.windowSeen,"Every frozen window binds a unique enrolled source");r.windowSeen=true;
			if(!r.eligible || mode.equals("support") && r.model!=supportModel){return;}require(Integer.parseInt(RrnaResourceIO.field(read.id,"start0"))==r.ts && Integer.parseInt(RrnaResourceIO.field(read.id,"stop0"))==r.te
				&& Integer.parseInt(RrnaResourceIO.field(read.id,"host"))==r.host && RrnaResourceIO.field(read.id,"split").equals(SPLITS[r.split]) && r.e<read.length() && r.te<read.length(),
				"Realign the original unshifted genomic window with the original host/truth coordinates");r.bases=read.bases;
		});for(Row r:rows){require(!r.eligible || mode.equals("support") && r.model!=supportModel || r.bases!=null,"Every included eligible member must retain its frozen window: "+r.id);}
	}
	private final class Worker extends Thread{
		Worker(AtomicInteger n){next=n;}
		@Override public void run(){try{final AlignmentStats stats=new AlignmentStats(true);final int[] pos=new int[2];
			for(int batch;(batch=next.getAndAdd(4))<rows.size();){for(int i=batch;i<Math.min(batch+4,rows.size());i++){
				final Row r=rows.get(i);if(!r.eligible || mode.equals("support") && r.model!=supportModel){continue;}final byte[][] v=variants[r.model];
				final float baseline=MitoEndpointAssignments.alignIndependent(library[r.model],r.bases,pos);
				require(pos[0]==r.s && pos[1]==r.e && baseline==r.identity,"Ordinary native baseline must reproduce the saved nontrace assignment exactly: "+r.id);
				r.result=new Fit[v==null?1:v.length];
				for(int j=0;j<r.result.length;j++){
					final byte[] model=v==null?library[r.model]:v[j];final Fit f=new Fit();r.result[j]=f;
					require(model.length<=r.bases.length,"Edited query exceeds the frozen genomic window; expand context explicitly before testing, never swap or drop: "+r.id+"/"+(v==null?"UNEDITED":variantNames[r.model][j]));
					if(v==null){f.identity=baseline;f.s=pos[0];f.e=pos[1];}
					else{stats.clear();stats.doTrace=true;f.identity=QuantumAligner.alignAndTraceStatic(model,r.bases,stats);f.s=stats.rStart;f.e=stats.rStop;
						require(stats.matchString!=null && stats.matchString.length>0,"Selected models require actual terminal traces");f.score=stats.score;f.ins=stats.ins;f.subs=stats.subs;
						final byte[] trace=stats.matchString;f.last=(char)trace[trace.length-1];f.first=(char)trace[0];f.trailI=0;f.leadI=0;
						for(int k=trace.length-1;k>=0 && trace[k]=='I';k--){f.trailI++;}for(int k=0;k<trace.length && trace[k]=='I';k++){f.leadI++;}
						f.tail=new String(trace,Math.max(0,trace.length-16),Math.min(16,trace.length),java.nio.charset.StandardCharsets.US_ASCII);
						f.repaired=stats.opposingIndelsWereRepaired();
						f.traceIdentity=stats.identity;f.ns=stats.ns;
						if(model[model.length-1]=='N' && !f.repaired && (f.last=='N' || f.last=='I')){
							final int prefixStop=f.last=='N'?f.e-1:f.e;
							if(prefixStop>=0 && prefixStop+1<r.bases.length){f.nMargin=terminalNMargin(model.length,r.bases.length,prefixStop);}
						}
						f.tieWitness=f.e==r.te-1 && f.trailI>0 && !f.repaired && terminalTie(model,r.bases,f.e);
					}
					require(f.s>=0 && f.e>=f.s && f.e<r.bases.length && Float.isFinite(f.identity),"Each alternative must return a valid inclusive full-model fit: "+r.id);
					if(j==0){require(f.s==r.s && f.e==r.e,"Unedited trace must preserve frozen coordinates; normalized trace identity is separate: "+r.id);}
				}
			}}success=true;
		}catch(Throwable t){failure=t;next.set(rows.size());}}
		final AtomicInteger next;boolean success;Throwable failure;
	}
	/** A valid equal-score alternate last operation: replace final I by S and
	 * consume the next genomic base. Excludes N and unequal final-column penalties.
	 * This witnesses a tied endpoint, not that the alternative was explored. */
	public static boolean terminalTie(byte[] q,byte[] ref,int stop){
		require(q!=null && ref!=null && q.length>0 && stop>=0 && stop<ref.length,"Tie evidence uses valid full-query inclusive bounds");
		if(stop+1>=ref.length){return false;}final byte a=q[q.length-1],b=ref[stop+1];
		if(a=='N' || b=='N' || a==b){return false;}
		// Mirrors QuantumAligner insPad and final-row insScoreMod, source lines171/255/280.
		final int pad=-16-ref.length/128-(int)Math.sqrt(ref.length)-(5*Math.max(0,ref.length-q.length))/4;
		return Math.max(0,(stop+1+pad-q.length)/2)==Math.max(0,(stop+2+pad-q.length)/2);
	}
	/** Same-prefix final N diagonal (0) versus final insertion (-1), including
	 * QuantumAligner:455 final-row positional penalties. Positive favors consuming
	 * the next reference base. This local witness does not prove global uniqueness. */
	public static int terminalNMargin(int qLen,int rLen,int prefixStop){
		require(qLen>0 && rLen>=qLen && prefixStop>=0 && prefixStop+1<rLen,"Compare valid adjacent endpoints for one final N on the same aligned prefix");
		final int pad=-16-rLen/128-(int)Math.sqrt(rLen)-(5*Math.max(0,rLen-qLen))/4;
		return 1+Math.max(0,(prefixStop+1+pad-qLen)/2)-Math.max(0,(prefixStop+2+pad-qLen)/2);
	}
	/** Real bases beyond a completely aligned model tip are kept separate from an
	 * unaligned terminal model insertion, so annotation offsets cannot manufacture
	 * extension support. No new base or edited model is generated here. */
	private void writeSupport(){
		final ObjectMap<String,long[]> groups=new ObjectMap<String,long[]>(String.class,long[].class);
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".support_members.tsv.gz"),summary=RrnaResourceIO.writer(prefix+".support_summary.tsv");
		final ByteBuilder b=new ByteBuilder();int n=0;
		try{out.println("locus\thost\tsplit\tinput_phylum\tmodel_index\tend\traw_error5\traw_error3\tmodel_length\tannotated_span\tnative_identity\ttrace_return_identity\tunpenalized_trace_score\tfirst_op\tlast_op\tleading_I\ttrailing_I\ttrace_tail\ttrace_normalized\tterminal_I_to_S_tie_witness\tstatus\tsupported_next_base\tannotated_overhang\tmodel_terminal12\tgene_terminal12");
			for(Row r:rows){if(!r.eligible || r.model!=supportModel){continue;}require(r.result!=null && r.result.length==1,"Support records must contain only one unedited trace per declared source");n++;
				final Fit f=r.result[0];final byte[] model=library[r.model];final boolean stop=supportEnd==3;
				final int distance=stop?r.te-f.e:f.s-r.ts;final char op=stop?f.last:f.first;final int inserted=stop?f.trailI:f.leadI;
				final int next=stop?f.e+1:f.s-1;final byte nextBase=distance>0?r.bases[next]:(byte)'N';
				final String status=supportStatus(f.repaired,inserted,op,stop && f.tieWitness,distance,model[stop?model.length-1:0],r.bases[stop?f.e:f.s],nextBase);
				char base='.';String overhang="NA";
				if(status.equals("ALIGNED_MODEL_WITH_ANNOTATED_OVERHANG") || status.equals("AMBIGUOUS_TIP_OR_OVERHANG")){
					require(next>=r.ts && next<=r.te,"Overhang evidence must be inside the original annotated gene");
					base=(char)nextBase;final int take=Math.min(12,distance);overhang=new String(r.bases,stop?next:next-take+1,take,java.nio.charset.StandardCharsets.US_ASCII);
				}
				final int error=stop?r.e-r.te:r.s-r.ts;
				for(int pass=0;pass<2;pass++){final String key=SPLITS[pass==0?0:r.split]+'\t'+status+'\t'+error+'\t'+base;long[] c=groups.get(key);if(c==null){c=new long[1];groups.put(key,c);}c[0]++;}
				final int mlen=Math.min(12,model.length),glen=Math.min(12,r.te-r.ts+1);
				final String mt=new String(model,stop?model.length-mlen:0,mlen,java.nio.charset.StandardCharsets.US_ASCII),gt=new String(r.bases,stop?r.te-glen+1:r.ts,glen,java.nio.charset.StandardCharsets.US_ASCII);
				out.print(b.clear().append(r.id).tab().append(r.host).tab().append(SPLITS[r.split]).tab().append(r.phylum).tab().append(r.model).tab().append(supportEnd)
					.tab().append(r.s-r.ts).tab().append(r.e-r.te).tab().append(model.length).tab().append(r.te-r.ts+1).tab().appendSlow(r.identity).tab().appendSlow(f.identity).tab().append(f.score)
					.tab().append(f.first).tab().append(f.last).tab().append(f.leadI).tab().append(f.trailI).tab().append(f.tail).tab().append(f.repaired).tab().append(f.tieWitness)
					.tab().append(status).tab().append(base).tab().append(overhang).tab().append(mt).tab().append(gt).nl());
			}
			require(n==expectedSupport,"Support analysis must preserve the exact previously counted model membership: "+n+" != "+expectedSupport);
			summary.println("split\tstatus\tsigned_error_at_selected_end\tsupported_next_base\tmembers");long all=0,split=0;
			for(String key:keys(groups)){final long c=groups.get(key)[0];if(key.startsWith("ALL\t")){all+=c;}else{split+=c;}summary.print(b.clear().append(key).tab().append(c).nl());}
			require(all==n && split==n,"Source-support statuses and host splits must independently conserve every member");
		}finally{RrnaResourceIO.finish(out,summary);}
		System.err.println("MITO_MODEL_SOURCE_SUPPORT_PASS model="+supportModel+" end="+supportEnd+" members="+n+" alignment_calls="+(2L*n)+" native_baseline_parity=true model_edits=0");
	}
	/** Canonical bases on both sides of the mapped model tip and the real next
	 * gene base are necessary evidence; Quantum's m/S alphabet alone is not. */
	public static String supportStatus(boolean repaired,int inserted,char op,boolean tie,int distance,byte modelTip,byte alignedTip,byte next){
		require(inserted>=0 && (op=='m' || op=='S' || op=='N' || op=='I' || op=='D'),"Support classification requires an actual native terminal operation");
		if(repaired){return "NORMALIZED_TRACE_UNRESOLVED";}
		if(inserted>0){return tie?"TERMINAL_I_EQUAL_SCORE_WITNESS":"UNMAPPED_MODEL_TIP";}
		if(op=='D'){return "TERMINAL_DELETION_UNRESOLVED";}
		if(distance==0){return "AT_ANNOTATED_END";}
		if(distance<0){return "PAST_ANNOTATED_END";}
		return op=='N' || RrnaResourceIO.base(modelTip)<0 || RrnaResourceIO.base(alignedTip)<0 || RrnaResourceIO.base(next)<0
			? "AMBIGUOUS_TIP_OR_OVERHANG":"ALIGNED_MODEL_WITH_ANNOTATED_OVERHANG";
	}
	private void writeTests(){
		RrnaResourceIO.fresh(prefix,".test_members.tsv.gz",".test_summary.tsv",".test_bases.tsv",".control.tsv");
		final ByteStreamWriter detail=RrnaResourceIO.writer(prefix+".test_members.tsv.gz"),summary=RrnaResourceIO.writer(prefix+".test_summary.tsv"),baseReport=RrnaResourceIO.writer(prefix+".test_bases.tsv"),control=RrnaResourceIO.writer(prefix+".control.tsv");
		final ObjectMap<String,Changes> bases=new ObjectMap<String,Changes>(String.class,Changes.class);
		final ByteBuilder b=new ByteBuilder();final Changes[][][] sums=new Changes[158][][];long unchanged=0,selectedN=0,alignments=0,alternativeRows=0;
		for(int m=5;m<158;m++){if(variants[m]!=null){sums[m]=new Changes[variants[m].length][3];}}
		try{detail.println("locus\thost\tsplit\tinput_phylum\tmodel_index\tmodel\tedited_end\tvariant\tmodel_length\ttruth_start0\ttruth_stop0\traw_error5\traw_error3\tnew_error5\tnew_error3\ttrace_return_identity\tunpenalized_trace_score\tinsertions\tsubstitutions\tfirst_op\tlast_op\tleading_I\ttrailing_I\ttrace_tail\ttrace_normalized\tterminal_I_to_S_tie_witness\tquery_terminal\ttruth_terminal\ttrace_stats_identity\ttrace_N_count\tterminal_N_consume_vs_insert_margin\tliteral_comparator_variant");
			for(Row r:rows){if(!r.eligible){continue;}alignments+=variants[r.model]==null?1:variants[r.model].length+1;if(variants[r.model]==null){unchanged++;continue;}selectedN++;alternativeRows+=r.result.length;
				for(int j=0;j<r.result.length;j++){final Fit f=r.result[j];final int a=f.s-r.ts,z=f.e-r.te;
					for(int split:new int[]{0,r.split}){Changes c=sums[r.model][j][split];if(c==null){c=new Changes();sums[r.model][j][split]=c;}c.add(r.s-r.ts,r.e-r.te,a,z);c.trace(f);}
					for(int pass=0;pass<2;pass++){final int split=pass==0?0:r.split;
						for(int group=0;group<BASE_GROUPS.length;group++){
							final Fit literal=r.result[literalComparator[r.model]];
							if((group==1 || group==2) && r.e-r.te!=-1 || group==2 && r.s!=r.ts
								|| group==3 && r.e!=r.te || group>=4 && (r.s!=r.ts || r.e!=r.te)
								|| group==5 && literal.s==r.ts && literal.e==r.te){continue;}
							final String key=r.model+"\t"+names[r.model]+'\t'+ends[r.model]+'\t'+variantNames[r.model][j]+'\t'+SPLITS[split]+'\t'+BASE_GROUPS[group]+'\t'+(char)r.bases[r.te];
							Changes c=bases.get(key);if(c==null){c=new Changes();bases.put(key,c);}c.add(r.s-r.ts,r.e-r.te,a,z);c.trace(f);
							if(z==0){c.correct3Ops[opIndex(f.last)]++;}
						}
					}
					b.clear().append(r.id).tab().append(r.host).tab().append(SPLITS[r.split]).tab().append(r.phylum).tab().append(r.model).tab().append(names[r.model])
						.tab().append(ends[r.model]).tab().append(variantNames[r.model][j]).tab().append(variants[r.model][j].length).tab().append(r.ts).tab().append(r.te)
						.tab().append(r.s-r.ts).tab().append(r.e-r.te).tab().append(a).tab().append(z).tab().appendSlow(f.identity).tab().append(f.score).tab().append(f.ins).tab().append(f.subs)
						.tab().append(f.first).tab().append(f.last).tab().append(f.leadI).tab().append(f.trailI).tab().append(f.tail).tab().append(f.repaired).tab().append(f.tieWitness)
						.tab().append((char)variants[r.model][j][ends[r.model]==5?0:variants[r.model][j].length-1]).tab().append((char)r.bases[ends[r.model]==5?r.ts:r.te])
						.tab().appendSlow(f.traceIdentity).tab().append(f.ns).tab();
					if(f.nMargin==Integer.MIN_VALUE){b.append("NA");}else{b.append(f.nMargin);}
					detail.print(b.tab().append(variantNames[r.model][literalComparator[r.model]]).nl());
				}
			}
			summary.println("model_index\tmodel\tedited_end\tvariant\tsplit\tmembers\traw_both_exact\tnew_both_exact\texact_gained\texact_lost\traw_within10\tnew_within10\twithin10_gained\twithin10_lost\tterminal_I_to_S_tie_witnesses\traw_exact5\tnew_exact5\traw_exact3\tnew_exact3\tabsolute5_improved\tabsolute5_worsened\tabsolute3_improved\tabsolute3_worsened\texact5_gained\texact5_lost\texact3_gained\texact3_lost\t"+NCOUNTS);
			for(int m=5;m<158;m++){if(variants[m]==null){continue;}for(int j=0;j<variants[m].length;j++){for(int split=0;split<3;split++){Changes c=sums[m][j][split];if(c==null){c=new Changes();}
				b.clear().append(m).tab().append(names[m]).tab().append(ends[m]).tab().append(variantNames[m][j]).tab().append(SPLITS[split]).tab().append(c.n)
					.tab().append(c.rawExact).tab().append(c.newExact).tab().append(c.gainExact).tab().append(c.lossExact).tab().append(c.rawNear).tab().append(c.newNear).tab().append(c.gainNear).tab().append(c.lossNear).tab().append(c.ties)
					.tab().append(c.raw5).tab().append(c.new5).tab().append(c.raw3).tab().append(c.new3).tab().append(c.better5).tab().append(c.worse5).tab().append(c.better3).tab().append(c.worse3)
					.tab().append(c.gain5).tab().append(c.loss5).tab().append(c.gain3).tab().append(c.loss3);c.appendN(b);summary.print(b.nl());}}}
			require(unchanged+selectedN==15084 && alignments==15084+alternativeRows,"Selected alternatives plus native baseline controls must conserve all eligible assignments");
			require(expectedSupport==0 || selectedN==expectedSupport,"Explicit selected-member denominator must match the completed model census");
			baseReport.println("model_index\tmodel\tedited_end\tvariant\tsplit\tbaseline_group\ttruth_3prime_base\tmembers\traw_both_exact\tnew_both_exact\texact_gained\texact_lost\traw_exact3\tnew_exact3\tcorrect3_terminal_m\tcorrect3_terminal_S\tcorrect3_terminal_I\tcorrect3_terminal_N\tcorrect3_terminal_D\traw_within10\tnew_within10\twithin10_gained\twithin10_lost\t"+NCOUNTS);
			for(String key:keys(bases)){final Changes c=bases.get(key);long operations=0;for(long n:c.correct3Ops){operations+=n;}require(operations==c.new3,"Every correct3'end must have exactly one observed terminal trace operation");
				b.clear().append(key).tab().append(c.n).tab().append(c.rawExact).tab().append(c.newExact).tab().append(c.gainExact).tab().append(c.lossExact).tab().append(c.raw3).tab().append(c.new3);
				for(long n:c.correct3Ops){b.tab().append(n);}b.tab().append(c.rawNear).tab().append(c.newNear).tab().append(c.gainNear).tab().append(c.lossNear);c.appendN(b);baseReport.print(b.nl());}
			control.println("eligible\tselected_models\tselected_members\tunchanged_model_members\tbaseline_coordinate_identity_mismatches\tactual_alignment_calls\tclaim");
			control.print(b.clear().append(15084).tab().append(selected).tab().append(selectedN).tab().append(unchanged).append("\t0\t").append(alignments).append("\tFROZEN_ASSIGNMENTS_NOT_MODEL_COMPETITION_OR_CALLER").nl());
		}finally{RrnaResourceIO.finish(detail,summary,baseReport,control);}
		System.err.println("MITO_MODEL_TERMINAL_TEST_PASS eligible=15084 selected_models="+selected+" selected_members="+selectedN+" unchanged_model_members="+unchanged+" alignments="+alignments+" baseline_parity=true");
	}
	private static boolean near(int a,int z,int t){return Math.abs((long)a)<=t && Math.abs((long)z)<=t;}
	private static int opIndex(char op){final int i="mSIND".indexOf(op);require(i>=0,"A trace must have a known final operation before cross-base attribution");return i;}
	private static <T> ArrayList<String> keys(ObjectMap<String,T> map){final ArrayList<String> a=new ArrayList<String>();for(String k:map.keys()){if(k!=null){a.add(k);}}Collections.sort(a);return a;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Counts{
		void add(int a,int z){n++;errors[0].increment(a);errors[1].increment(z);if(a==0){exact5++;}if(z==0){exact3++;}if(a==0 && z==0){exact++;}if(near(a,z,10)){within10++;}
			if(a==0 && z==-1){short3++;}if(a==0 && z==1){long3++;}if(a==-1 && z==0){long5++;}if(a==1 && z==0){short5++;}
			if(a==-1){m5minus++;}if(a==1){m5plus++;}if(z==-1){m3minus++;}if(z==1){m3plus++;}
		}
		void append(ByteBuilder b){b.tab().append(n).tab().append(exact5).tab().append(exact3).tab().append(exact).tab().append(within10).tab().append(short3).tab().append(long3).tab().append(long5).tab().append(short5).tab().append(m5minus).tab().append(m5plus).tab().append(m3minus).tab().append(m3plus);}
		long n,exact5,exact3,exact,within10,short3,long3,long5,short5,m5minus,m5plus,m3minus,m3plus;
		final IntHashMap[] errors={new IntHashMap(8),new IntHashMap(8)};
	}
	private static final class Changes{
		void trace(Fit f){if(f.tieWitness){ties++;}if(f.nMargin==Integer.MIN_VALUE){return;}
			if(f.last=='N'){nConsumed++;}else{require(f.last=='I',"N local-margin evidence requires a final N or I");nInserted++;}
			if(f.nMargin>0){nPositive++;}else if(f.nMargin==0){nTied++;}else{nNegative++;}}
		void appendN(ByteBuilder b){require(nConsumed+nInserted==nPositive+nTied+nNegative,"Every local N comparison has exactly one operation and score-margin sign");
			b.tab().append(nConsumed).tab().append(nInserted).tab().append(nPositive).tab().append(nTied).tab().append(nNegative);}
		void add(int a,int z,int x,int y){n++;final boolean r=a==0 && z==0,s=x==0 && y==0,u=near(a,z,10),v=near(x,y,10);
			if(a==0){raw5++;}if(x==0){new5++;}if(z==0){raw3++;}if(y==0){new3++;}
			if(a!=0 && x==0){gain5++;}if(a==0 && x!=0){loss5++;}if(z!=0 && y==0){gain3++;}if(z==0 && y!=0){loss3++;}
			if(Math.abs((long)x)<Math.abs((long)a)){better5++;}if(Math.abs((long)x)>Math.abs((long)a)){worse5++;}
			if(Math.abs((long)y)<Math.abs((long)z)){better3++;}if(Math.abs((long)y)>Math.abs((long)z)){worse3++;}
			if(r){rawExact++;}if(s){newExact++;}if(!r && s){gainExact++;}if(r && !s){lossExact++;}if(u){rawNear++;}if(v){newNear++;}if(!u && v){gainNear++;}if(u && !v){lossNear++;}}
		long n,rawExact,newExact,gainExact,lossExact,rawNear,newNear,gainNear,lossNear,ties,raw5,new5,raw3,new3,better5,worse5,better3,worse3,gain5,loss5,gain3,loss3;
		long nConsumed,nInserted,nPositive,nTied,nNegative;
		final long[] correct3Ops=new long[5];
	}
	private static final class Row{String id,phylum,disposition;int host,split,model,ts,te,s,e;float identity;boolean seen,eligible,windowSeen;byte[] bases;Fit[] result;}
	private static final class Fit{int s,e,score,ins,subs,leadI,trailI,ns,nMargin=Integer.MIN_VALUE;float identity,traceIdentity;char first,last;String tail;boolean repaired,tieWitness;}
	private String mode="counts",assignments,enrollment,consensus,corpus,edits,prefix;private byte[][] library;private String[] names;private byte[][][] variants;private final int[] ends=new int[158];private int selected,supportModel=-1,supportEnd=3,expectedSupport;
	private final ObjectMap<String,Row> byId=new ObjectMap<String,Row>(String.class,Row.class);private final ArrayList<Row> rows=new ArrayList<Row>();
	private final String[][] variantNames=new String[158][];
	private final int[] literalComparator=new int[158];
	private static final byte[] BASES={'A','C','G','T'};
	private static final String[] SPLITS={"ALL","TRAIN","VALIDATION"};
	private static final String[] BASE_GROUPS={"ALL","RAW_3_MINUS1","RAW_0_MINUS1","RAW_3_EXACT","RAW_BOTH_EXACT","LITERAL_BOTH_EXACT_LOSS"};
	private static final String NCOUNTS="local_witness_N_consumed\tlocal_witness_N_inserted\tlocal_witness_consume_margin_positive\tlocal_witness_consume_margin_zero\tlocal_witness_consume_margin_negative";
	private static final String COUNTS="members\traw_exact5\traw_exact3\traw_both_exact\traw_both_within10\tpair_0_minus1\tpair_0_plus1\tpair_minus1_0\tpair_plus1_0\tmarginal5_minus1\tmarginal5_plus1\tmarginal3_minus1\tmarginal3_plus1";
	private static final String HEADER="locus_id\tarm\thost\tsplit\tsource_status\tpilot_caller_status\tdisposition\tmodel_index\tmodel\tshared9\tmodels_aligned\tassignment_identity_shorter_query\tgenomic_identity_full_model\ttruth_start0\ttruth_stop0\traw_start0\traw_stop0\tcorrection5\tcorrection3\tpartial\tpseudo";
}
