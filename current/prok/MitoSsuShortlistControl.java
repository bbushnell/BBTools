package prok;

import java.util.ArrayList;
import java.util.Collections;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Bounded native16S topN control. Reproduces saved model/window rows before
 * interpreting an expanded shortlist. Isolated family candidates, not final DP calls.
 * @author Brian Bushnell, Raiden */
public final class MitoSsuShortlistControl {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, MitoSsuShortlistControl.class, false);
		String in=null, stages=null, contig=null, out=null, cases=null;
		int strand=-1, truthStart=-1, truthStop=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq), value=arg.substring(eq+1);
			if(key.equals("in")){in=value;}else if(key.equals("stages")){stages=value;}
			else if(key.equals("contig")){contig=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("strand")){strand=Integer.parseInt(value);}
			else if(key.equals("cases")){cases=value;}
			else if(key.equals("truthstart")){truthStart=Integer.parseInt(value);}
			else if(key.equals("truthstop")){truthStop=Integer.parseInt(value);}
			else{throw new IllegalArgumentException(key);}
		}
		require(stages!=null && out!=null && (cases==null ? in!=null && contig!=null && (strand==0 || strand==1)
			: in==null && contig==null && strand==-1),
			"Bind one source contig/strand or the explicit22-case evidence table, plus saved stages and fresh output");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.seedMinHits=2;
		final NcrnaFamily family=NcrnaFamilyLoader.load("16S", options, new FamilyResourceResolver());
		NcrnaFamilyRegistry.enablePhysicalEdges(family);//CallGenes default ncrnaedge phase precedes worker construction.
		require(MitoSsuModelPolicy.expanded(family.name, family.modelNames) && family.indexK==9
			&& family.indexTopN==1 && !family.adaptive && family.fixedMinHits==1,
			"Control must use the expanded ten-model16S resource and actual K9/topN1 configuration");
		if(cases!=null){
			require(truthStart==-1 && truthStop==-1, "Batch truth comes only from the bound case table");
			batch(family, cases, stages, out);Shared.closeStream(pp.outstream);return;
		}
		final String wanted=contig;final Read[] selected=new Read[1];
		RrnaResourceIO.stream(in, r->{if(RrnaResourceIO.first(r.id).equals(wanted)){
			require(selected[0]==null && r.bases!=null, "Exactly one original contig is required");
			selected[0]=r;
		}});
		require(selected[0]!=null, "Requested contig absent from original FASTA");
		if(strand==1){AminoAcid.reverseComplementBasesInPlace(selected[0].bases);}
		ArrayList<Case> exact=null;
		if(truthStart!=-1 || truthStop!=-1){
			require(truthStart>0 && truthStop>=truthStart && truthStop<=selected[0].length(), "Explicit single-control truth must fit its source");
			exact=new ArrayList<Case>();final Case c=new Case("prok_control",contig,in,strand,truthStart,truthStop,family.library.length);
			c.bases=selected[0].bases;exact.add(c);
		}
		for(int top:new int[]{1, family.library.length}){
			run(family, selected[0], strand, top, out+".top"+top, false, top==1?null:exact);
			if(top==1){
				final ArrayList<String> before=rows(stages, contig, strand);
				final ArrayList<String> after=rows(out+".top1.stages.tsv", contig, strand);
				require(!before.isEmpty() && before.equals(after),
					"TopN1 must reproduce every actual model/window row before interpreting the control");
			}
		}
		if(exact!=null){writeModels(family,exact,out);}
		System.err.println("MITO_SSU_SHORTLIST_CONTROL_PASS production_top1_model_windows_exact isolated_family_no_DP");
		Shared.closeStream(pp.outstream);
	}

	static void run(NcrnaFamily f, Read r, int strand, int top, String prefix, boolean oracle, ArrayList<Case> cases){
		assert(top>=1 && top<=f.library.length):"Only the observed one-model policy and bounded full shortlist are compared";
		RrnaResourceIO.fresh(prefix, ".stages.tsv", ".attempts.tsv", ".candidates.gff");
		final NcrnaScavenger scav=new NcrnaScavenger(f.library, f.models, f.modelNames,
			f.kmerSet, f.kLong, f.minLen, f.windowPad, f.indexK, top, f.adaptive,
			f.adaptFloor, f.adaptTopFrac, f.adaptQFrac, f.fixedMinHits,
			f.scoreA, f.scoreB, f.idPass, f.idBorderline);
		GeneCaller.applyNcrnaFamilyControls(scav, f);
		scav.setMappingOracleExhaustive(oracle);
		final ByteStreamWriter stages=RrnaResourceIO.writer(prefix+".stages.tsv");
		final ByteStreamWriter attempts=RrnaResourceIO.writer(prefix+".attempts.tsv");
		final ByteStreamWriter calls=RrnaResourceIO.writer(prefix+".candidates.gff");
		final ByteBuilder b=new ByteBuilder();
		stages.println(NcrnaStageDiagLogSink.HEADER);
		attempts.println("contig\tstrand\tpass\twindowStart0\twindowStop0\trank\tmodel\tsharedK9\tidentity\talignedLength\tidentityAndSpanPass");
		scav.setStageDiagSink(new NcrnaStageDiagLogSink(stages), "16S");
		scav.setModelAttemptSink((name, s, pass, start, stop, rank, model, shared, id, length, accepted)->{
			assert(model>=0 && model<f.modelNames.length):"Attempt rows bind actual native model indices";
			attempts.print(b.clear().append(RrnaResourceIO.first(name)).tab().append(s).tab().append(pass)
				.tab().append(start).tab().append(stop).tab().append(rank).tab().append(f.modelNames[model])
				.tab().append(shared).tab().appendSlow(id).tab().append(length).tab().append(accepted).nl());
			if(cases!=null){for(Case c:cases){
				final int a=s==0 ? c.start-1 : r.length()-c.stop;
				final int z=s==0 ? c.stop-1 : r.length()-c.start;
				if(start>a || stop<z){continue;}//Only actual windows containing the entire truth interval.
				c.attempts[model]++;c.eligible[model]|=accepted;
				if(id>c.identity[model]){c.identity[model]=id;c.length[model]=length;c.ws[model]=start;c.we[model]=stop;}
			}}
		});
		try{
			final ArrayList<Orf> result=scav.scavenge(r.id, r.bases, strand, new ArrayList<int[]>());
			calls.println("##gff-version 3");
			for(Orf orf:result){if(strand==1){orf.flip();}calls.println(orf.toGff());}
		}finally{RrnaResourceIO.finish(stages, attempts, calls);}
	}

	/** Compare native all-model verification on the18existing genomes, conserving22truth loci.
	 * The second policy bypasses only mapping/early exit; seed/window/verification/ends stay native. */
	static void batch(NcrnaFamily family, String casesPath, String stages, String prefix){
		final ArrayList<Case> all=new ArrayList<Case>();
		final map.ObjectMap<String,ArrayList<Case>> groups=new map.ObjectMap(String.class, ArrayList.class);
		final ByteFile in=RrnaResourceIO.open(casesPath);final LineParser1 p=new LineParser1('\t');
		try{
			p.set(in.nextLine());require(p.terms()>=31 && p.termEquals("kind",1) && p.termEquals("source_fasta",26), "Expected source-bound evidence case schema");
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()>=31, "Evidence case row is incomplete");
				if(!p.termEquals("ORGANELLE_SSU_TRUTH",1)){continue;}
				final Case c=new Case(p, family.library.length);all.add(c);
				ArrayList<Case> group=groups.get(c.contig);
				if(group==null){group=new ArrayList<Case>();groups.put(c.contig,group);}
				else{require(group.get(0).fasta.equals(c.fasta) && group.get(0).strand==c.strand,
					"This bounded replay requires one source/true strand per organelle contig");}
				group.add(c);
			}
		}finally{require(!in.close(), "Evidence case table read failed");}
		require(all.size()==22 && groups.size()==18, "Preserve all22truth loci on exactly18original organelle contigs");
		final ArrayList<String> order=new ArrayList<String>();for(String name:groups.keys()){if(name!=null){order.add(name);}}Collections.sort(order);
		for(String contig:order){
			final ArrayList<Case> cases=groups.get(contig);final Case first=cases.get(0);final Read[] read=new Read[1];
			RrnaResourceIO.stream(first.fasta,r->{if(RrnaResourceIO.first(r.id).equals(contig)){
				require(read[0]==null && r.bases!=null, "Original organelle contig must occur once");read[0]=r;
			}});
			require(read[0]!=null, "Missing original organelle contig "+contig);
			for(Case c:cases){require(c.stop<=read[0].length(), "Truth interval must fit original source sequence");}
			if(first.strand==1){AminoAcid.reverseComplementBasesInPlace(read[0].bases);}
			for(Case c:cases){c.bases=read[0].bases;}
			final String stem=prefix+"."+contig;
			run(family,read[0],first.strand,1,stem+".top1",false,null);
			final ArrayList<String> expected=rows(stages,contig,first.strand);
			require(expected.equals(rows(stem+".top1.stages.tsv",contig,first.strand)),
				"Every top1 model/window stream must reproduce production before all-model interpretation: "+contig);
			run(family,read[0],first.strand,family.library.length,stem+".allmodels",true,cases);
		}
		writeModels(family,all,prefix);
		System.err.println("MITO_SSU_ALL_MODEL_CONTROL_PASS contigs=18 loci=22 reported_pairs=132 top1_parity=true isolated_family_no_DP");
	}
	static void writeModels(NcrnaFamily family,ArrayList<Case> all,String prefix){
		require(!idaligner.GlocalAligner.GLOBAL, "Exact control requires glocal, not whole-reference global mode");
		RrnaResourceIO.fresh(prefix,".models.tsv");final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".models.tsv");
		final ByteBuilder b=new ByteBuilder();int count=0;
		try{
			out.println("# Maximum identity across actual native attempts in windows containing the entire locus, including physical-edge retry; eligibility precedes endpoint/path selection. This is not exact-gene global identity or final-caller recall. All-model raw ranks are library traversal positions.");
			out.println("# Quantum/Glocal controls both align the complete model as query to the SAME saved window; Glocal also aligns model to exact truth gene. Non-null stats prevent query/reference swapping. Returned identity is kept separate from traceback identity. Clipped native retries can differ from the full-model Quantum control.");
			out.println("case_id\tcontig\tstart1\tstop1\tstrand\tmodel\tnativeWindowIdentityMax\tnativeAlignedLengthAtMaximum\twindowStart0\twindowStop0\tcontainingWindowAttempts\tanyNativeWindowEligible\tcomparisonStatus\tquantumWindow\tglocalWindow\tglocalExactGene\tquantumTraceIdentity\tglocalWindowTraceIdentity\tglocalGeneTraceIdentity\tquantumWindowStart0\tquantumWindowStop0\tglocalWindowStart0\tglocalWindowStop0\tquantumTruthOverlapBp\tglocalTruthOverlapBp\ttruthLength\tmodelLength\tquantumTruthCoverage\tglocalTruthCoverage\trowRole");
			for(Case c:all){for(int m=0;m<family.library.length;m++){
				final boolean prokControl=c.id.equals("prok_control"),bacteria=prokControl && m==1;
				if(m!=0 && m<5 && !bacteria){continue;}//Six requested models plus the known bacterial positive control.
				final int a=c.strand==0?c.start-1:c.bases.length-c.stop, z=c.strand==0?c.stop:c.bases.length-c.start+1;
				final byte[] gene=java.util.Arrays.copyOfRange(c.bases,a,z);
				final idaligner.AlignmentStats x=new idaligner.AlignmentStats(true);
				final float xi=idaligner.GlocalAligner.alignAndTraceStatic(family.library[m],gene,x);
				require(Float.isFinite(xi), "Exact-gene comparison must run even when production had no containing window");
				b.clear().append(c.id).tab().append(c.contig).tab().append(c.start).tab().append(c.stop)
					.tab().append(c.strand==0?'+':'-').tab().append(family.modelNames[m]);
				if(c.attempts[m]==0){
					b.append("\tNA\tNA\tNA\tNA\t0\tNA\tNOT_ATTEMPTED_IN_CONTAINING_WINDOW\tNA\tNA\t").appendSlow(xi)
						.append("\tNA\tNA\t").appendSlow(x.identity).append("\tNA\tNA\tNA\tNA\tNA\tNA\t")
						.append(gene.length).tab().append(family.library[m].length).append("\tNA\tNA");
				}else{
					final byte[] window=java.util.Arrays.copyOfRange(c.bases,c.ws[m],c.we[m]+1);
					final idaligner.AlignmentStats q=new idaligner.AlignmentStats(true),g=new idaligner.AlignmentStats(true);
					final float qi=idaligner.QuantumAligner.alignAndTraceStatic(family.library[m],window,q);
					final float gi=idaligner.GlocalAligner.alignAndTraceStatic(family.library[m],window,g);
					require(Float.isFinite(qi) && Float.isFinite(gi), "Same-window comparisons must complete with finite native identities");
					if(c.ws[m]>0 && c.we[m]<c.bases.length-1){require(Float.floatToIntBits(qi)==Float.floatToIntBits(c.identity[m]), "Interior native Quantum maximum must reproduce on the exact model/window pair");}
					final int qa=c.ws[m]+q.rStart,qz=c.ws[m]+q.rStop,ga=c.ws[m]+g.rStart,gz=c.ws[m]+g.rStop;
					final int qo=overlap(qa,qz,a,z-1),go=overlap(ga,gz,a,z-1);
					b.tab().appendSlow(c.identity[m]).tab().append(c.length[m]).tab().append(c.ws[m]).tab().append(c.we[m])
						.tab().append(c.attempts[m]).tab().append(c.eligible[m]).append("\tSAME_MODEL_WINDOW\t")
						.appendSlow(qi).tab().appendSlow(gi).tab().appendSlow(xi).tab().appendSlow(q.identity)
						.tab().appendSlow(g.identity).tab().appendSlow(x.identity).tab().append(qa).tab().append(qz)
						.tab().append(ga).tab().append(gz).tab().append(qo).tab().append(go).tab().append(gene.length)
						.tab().append(family.library[m].length).tab().append(qo/(double)gene.length,9).tab().append(go/(double)gene.length,9);
				}
				out.print(b.tab().append(bacteria?"BACTERIAL_POSITIVE_CONTROL":prokControl?"PROK_CONTROL":"REQUESTED_ORGANELLE_PAIR").nl());count++;
			}}
			final int positiveControls=all.size()==1 && all.get(0).id.equals("prok_control")?1:0;
			require(count==6*all.size()+positiveControls, "Six requested models per locus plus the explicitly tagged bacterial positive control must be retained");
		}finally{RrnaResourceIO.finish(out);}
	}
	static int overlap(int a,int z,int start,int stop){return z<a?0:Math.max(0,Math.min(z,stop)-Math.max(a,start)+1);}
	static final class Case{
		Case(LineParser1 p,int models){
			this(p.parseString(0),p.parseString(3),p.parseString(26),strand(p.parseString(6)),p.parseInt(4),p.parseInt(5),models);
		}
		Case(String id_,String contig_,String fasta_,int strand_,int start_,int stop_,int models){
			id=id_;contig=contig_;fasta=fasta_;strand=strand_;start=start_;stop=stop_;
			require(start>0 && stop>=start && contig.indexOf('/')<0 && contig.indexOf('\\')<0, "Valid source geometry and safe accession token required");
			identity=new float[models];java.util.Arrays.fill(identity,-1f);length=new int[models];ws=new int[models];we=new int[models];attempts=new int[models];eligible=new boolean[models];
		}
		static int strand(String s){require(s.equals("+") || s.equals("-"), "Explicit true strand required");return s.equals("+")?0:1;}
		final String id,contig,fasta;final int start,stop,strand;final float[] identity;final int[] length,ws,we,attempts;final boolean[] eligible;
		byte[] bases;
	}

	static ArrayList<String> rows(String path, String contig, int strand){
		assert(path!=null && contig!=null):"Parity compares actual immutable stage streams";
		final ArrayList<String> result=new ArrayList<String>();
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, NcrnaStageDiagLogSink.HEADER, path);
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==21, "Native stage schema changed");
				if((p.termEquals("MODEL", 0) || p.termEquals("WINDOW", 0)) && p.termEquals("16S", 1)
					&& RrnaResourceIO.first(p.parseString(2)).equals(contig) && p.parseInt(3)==strand){
					result.add(new String(line, java.nio.charset.StandardCharsets.UTF_8));
				}
			}
		}finally{require(!in.close(), "Saved stage read failed");}
		Collections.sort(result);return result;
	}
	static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
}
