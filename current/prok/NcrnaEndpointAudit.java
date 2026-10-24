package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Describes signed endpoint offsets from one saved frozen-locus matching audit.
 * Positive offsets move toward the reference molecule's 3-prime direction.
 * Each stored offset is independently checked against its paired coordinates.
 * No calling, matching, truth changes or biological adjudication.
 * @author Raiden */
public final class NcrnaEndpointAudit {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaEndpointAudit.class,false);
		String in=null,out=null,arm=null,panel=null,compartment=null,family=null,geometry=null;int expected=-1,e5=-1,e3=-1,both=-1,selected=-1;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("in")){in=v;}else if(k.equals("outprefix")){out=v;}else if(k.equals("arm")){arm=v;}else if(k.equals("panel")){panel=v;}
			else if(k.equals("compartment")){compartment=v;}else if(k.equals("family")){family=v;}else if(k.equals("expected")){expected=Integer.parseInt(v);}
			else if(k.equals("exact5")){e5=Integer.parseInt(v);}else if(k.equals("exact3")){e3=Integer.parseInt(v);}else if(k.equals("both")){both=Integer.parseInt(v);}
			else if(k.equals("geometry")){geometry=v;}else if(k.equals("expectedselected")){selected=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown flag: "+k);}}
		require(in!=null && out!=null && arm!=null && panel!=null && compartment!=null && family!=null && expected>0 && e5>=0 && e3>=0 && both>=0,"Explicit selection and sealed count expectations required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		require(geometry==null ? selected==-1 : selected>0 && (geometry.equals("FULLSPAN_INTACT") || geometry.equals("EDGE") || geometry.equals("PARTIAL_INTERNAL") || geometry.equals("UNASSESSABLE")), "Geometry selection needs a recognized frozen stratum and explicit positive expectedselected count");
		final NcrnaEndpointAudit x=new NcrnaEndpointAudit();x.load(in,arm,panel,compartment,family,expected,e5,e3,both);
		if(geometry!=null){x.selectGeometry(geometry, selected);}x.write(out);Shared.closeStream(pp.outstream);
	}
	/** Validate the entire parent audit first: excluded pairs cannot hide corrupt offsets. */
	void selectGeometry(String geometry, int expected){
		assert(!rows.isEmpty()):"Geometry selection follows complete one-to-one audit validation";
		final int parent=rows.size();rows.removeIf(r->!r.geometry.equals(geometry));
		require(rows.size()==expected, "Frozen geometry count differs from the independently graded stratum: "+geometry+" observed="+rows.size()+" expected="+expected);
		System.err.println("ENDPOINT_GEOMETRY_SELECTED geometry="+geometry+" matched="+rows.size()+" parent_matched="+parent);
	}
	void load(String path,String arm,String panel,String compartment,String family,int expected,int exact5,int exact3,int exactBoth){
		assert(rows.isEmpty()):"One saved selection supplies every matched pair";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectMap<String,Pair> calls=new ObjectMap<String,Pair>(String.class,Pair.class);
		final ObjectSet<String> truthIds=new ObjectSet<String>(String.class);int ncall=0,n5=0,n3=0,nboth=0;
		try{RrnaResourceIO.header(in,HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==18,"Saved audit schema changed");
			if(!p.termEquals(arm,0) || !p.termEquals(panel,1) || !p.termEquals(compartment,4) || !p.termEquals(family,5) || !p.termEquals("TP",12)){continue;}
			require(p.termEquals("FROZEN_LOCUS_TRUTH_V1",17),"Only the approved frozen-locus matching audit is supported");
			final String assembly=p.parseString(2),contig=p.parseString(3),id=p.parseString(7),partner=p.parseString(13);
			if(p.termEquals("TRUTH",6)){
				final Pair r=new Pair(p);require(truthIds.add(assembly+'\t'+contig+'\t'+id) && calls.put(assembly+'\t'+contig+'\t'+partner,r)==null,"Duplicate reference or consumed call in one-to-one audit");
				rows.add(r);
			}else{
				require(p.termEquals("CALL",6),"Unknown matched audit role");final Pair r=calls.get(assembly+'\t'+contig+'\t'+id);
				require(r!=null && !r.complete && r.id.equals(partner) && p.termEquals(r.strand,11),"Matched CALL must follow its unique same-strand TRUTH row, as emitted by NcrnaClippingComparison.grade");
				r.callStart=p.parseInt(9);r.callEnd=p.parseInt(10);require(r.callStart>0 && r.callEnd>=r.callStart,"Paired caller coordinates must be nonempty genomic intervals");
				final boolean minus=r.strand.equals("-");final int d5=minus?r.end-r.callEnd:r.callStart-r.start,d3=minus?r.start-r.callStart:r.callEnd-r.end;
				require(d5==r.d5 && d3==r.d3 && p.parseInt(14)==d5 && p.parseInt(15)==d3,"Saved signed offsets must equal paired genomic coordinates in transcript orientation");
				r.complete=true;ncall++;if(d5==0){n5++;}if(d3==0){n3++;}if(d5==0 && d3==0){nboth++;}
			}
		}}finally{require(!in.close(),"Saved matching audit read failed");}
		require(rows.size()==expected && ncall==expected && n5==exact5 && n3==exact3 && nboth==exactBoth,"Selected matching/endpoint counts must equal the sealed grade summary");
		for(Pair r:rows){require(r.complete,"Missing paired CALL row");}System.err.println("ENDPOINT_AUDIT_COMPLETE matched="+rows.size()+" exact5="+n5+" exact3="+n3+" both="+nboth+" no_caller_runs=true");
	}
	void write(String prefix){
		assert(!rows.isEmpty()):"A validated paired selection is required before reporting offsets";
		RrnaResourceIO.fresh(prefix,".pairs.tsv",".histogram.tsv",".quantiles.tsv",".summary.tsv");
		final ByteStreamWriter pairs=RrnaResourceIO.writer(prefix+".pairs.tsv"),hist=RrnaResourceIO.writer(prefix+".histogram.tsv"),q=RrnaResourceIO.writer(prefix+".quantiles.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();final int[][] offsets=new int[2][rows.size()];final long[] states=new long[4];
		try{
			pairs.println("assembly\tcontig\tstrand\ttruth_id\tcall_id\ttruth_start\ttruth_end\tcall_start\tcall_end\td5\td3\tgeometry");
			for(int i=0;i<rows.size();i++){final Pair r=rows.get(i);offsets[0][i]=r.d5;offsets[1][i]=r.d3;states[(r.d5==0?1:0)+(r.d3==0?2:0)]++;
				pairs.print(b.clear().append(r.assembly).tab().append(r.contig).tab().append(r.strand).tab().append(r.id).tab().append(r.partner).tab()
					.append(r.start).tab().append(r.end).tab().append(r.callStart).tab().append(r.callEnd).tab().append(r.d5).tab().append(r.d3).tab().append(r.geometry).nl());}
			hist.println("end\tmeasure\toffset_bp\tloci");q.println("end\tmeasure\tpopulation\tloci\tminimum\tp10\tp50\tp90\tp95\tp99\tmaximum");
			for(int end=0;end<2;end++){final String label=end==0?"5prime":"3prime";
				for(int measure=0;measure<2;measure++){final String kind=measure==0?"signed":"absolute";final int[] values=offsets[end].clone();
					if(measure==1){for(int i=0;i<values.length;i++){require(values[i]!=Integer.MIN_VALUE,"Absolute offsets must fit an integer");values[i]=Math.abs(values[i]);}}
					Arrays.sort(values);for(int i=0;i<values.length;){int j=i+1;while(j<values.length && values[j]==values[i]){j++;}hist.print(b.clear().append(label).tab().append(kind).tab().append(values[i]).tab().append(j-i).nl());i=j;}
					quantiles(q,b,label,kind,"all_matched",values);
					int nonzero=0;for(int value:values){if(value!=0){nonzero++;}}final int[] errors=new int[nonzero];int i=0;for(int value:values){if(value!=0){errors[i++]=value;}}
					quantiles(q,b,label,kind,"nonexact_only",errors);
				}
			}
			summary.println("outcome\tloci\tmatched_total");final String[] names={"neither_exact","only_5prime_exact","only_3prime_exact","both_exact"};
			for(int i=0;i<states.length;i++){summary.print(b.clear().append(names[i]).tab().append(states[i]).tab().append(rows.size()).nl());}
		}finally{RrnaResourceIO.finish(pairs,hist,q,summary);}
	}
	static void quantiles(ByteStreamWriter out,ByteBuilder b,String end,String measure,String population,int[] a){
		assert(a!=null):"Quantiles use nearest rank on the sorted selected population";
		b.clear().append(end).tab().append(measure).tab().append(population).tab().append(a.length);
		if(a.length==0){b.append("\tNA\tNA\tNA\tNA\tNA\tNA\tNA");}
		else{b.tab().append(a[0]);for(int percentile:new int[]{10,50,90,95,99}){b.tab().append(a[(int)(((long)percentile*a.length+99)/100)-1]);}b.tab().append(a[a.length-1]);}out.print(b.nl());
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Pair{
		Pair(LineParser1 p){assembly=p.parseString(2);contig=p.parseString(3);id=p.parseString(7);partner=p.parseString(13);start=p.parseInt(9);end=p.parseInt(10);strand=p.parseString(11);d5=p.parseInt(14);d3=p.parseInt(15);geometry=p.parseString(16);
			require(start>0 && end>=start && (strand.equals("+") || strand.equals("-")) && !partner.equals("NA"),"Matched truth requires valid coordinates, strand and partner");}
		final String assembly,contig,id,partner,strand,geometry;final int start,end,d5,d3;int callStart,callEnd;boolean complete;
	}
	final ArrayList<Pair> rows=new ArrayList<Pair>();
	public static final String HEADER="arm\tpanel\tassembly\tcontig\tcompartment\tview\trole\tid\tlabel\tstart\tstop\tstrand\toutcome\tpartner\td5\td3\tgeometry\treference_kind";
}
