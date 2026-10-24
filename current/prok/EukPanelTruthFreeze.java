package prok;

import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Marks an independently accepted disposition table as explicitly versioned truth.
 * The committed producer must bind Ganyu's acceptance before invoking this tool.
 * Every data byte other than the single status field is preserved.
 * @author Raiden */
public final class EukPanelTruthFreeze {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,EukPanelTruthFreeze.class,false);String in=null,out=null,policy=EukPanelTruthDisposition.LEGACY_POLICY;long expected=-1;int version=1;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("in")){in=v;}else if(k.equals("out")){out=v;}else if(k.equals("expectedloci")){expected=Long.parseLong(v);}else if(k.equals("policy")){policy=v;}else if(k.equals("version")){version=Integer.parseInt(v);}else{throw new IllegalArgumentException("Unknown flag: "+k);}}
		require(in!=null && out!=null && expected>0,"Explicit reviewed input, fresh output and expected count required");Shared.setThreads(1);ReadWrite.setZipThreads(1);run(in,out,expected,policy,version);Shared.closeStream(pp.outstream);
	}
	static void run(String input,String output,long expected){
		run(input,output,expected,EukPanelTruthDisposition.LEGACY_POLICY);
	}
	static void run(String input,String output,long expected,String policy){
		run(input,output,expected,policy,1);
	}
	static void run(String input,String output,long expected,String policy,int version){
		assert(expected>0):"Freeze changes status only; the accepted locus count must be supplied";
		validateVersion(policy,version);
		final String candidate=version==1?"FREEZE_CANDIDATE_NOT_FROZEN":V2_CANDIDATE;
		RrnaResourceIO.fresh(output,"");final ByteFile in=RrnaResourceIO.open(input);final ByteStreamWriter out=RrnaResourceIO.writer(output);
		final LineParser1 p=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();final ObjectSet<String> seen=new ObjectSet<String>(String.class);long n=0;
		try{RrnaResourceIO.header(in,HEADER,input);out.println(HEADER);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==35 && p.termEquals(candidate,32) && p.termEquals(policy,34),"Only the accepted candidate of the explicitly selected policy/version can be frozen");
			require(seen.add(p.parseString(0)),"Duplicate locus ID cannot be frozen");
			p.setBounds(32);b.clear().append(line,0,p.a()).append("TRUTH_V").append(version).append(line,p.b(),line.length-p.b());out.print(b.nl());n++;
		}require(n==expected,"Freeze must preserve the reviewed locus denominator");}
		finally{final boolean error=in.close();RrnaResourceIO.finish(out);require(!error,"Frozen truth input read failed");}
		System.err.println("TRUTH_FREEZE_COMPLETE loci="+n+" status_only_change=true");
	}
	static void validateVersion(String policy,int version){
		require(version==1 && (policy.equals(EukPanelTruthDisposition.LEGACY_POLICY) || policy.equals(EukPanelTruthDisposition.PROK_POLICY))
			|| version==2 && policy.equals(V2_POLICY),"Unknown reviewed truth policy/version combination");
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String HEADER=EukPanelTruthDisposition.OUTPUT_HEADER;
	public static final String V2_POLICY="HISTORICAL_CLASS_ADDENDUM_v2",V2_CANDIDATE="FREEZE_CANDIDATE_V2_NOT_FROZEN";
}
