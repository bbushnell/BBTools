package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Aggregates saved final-member identities without repeating any alignment.
 * Means use the nine-decimal identities in the frozen member table. Partial and
 * pseudo strata may overlap; all source copies remain in the ALL denominator.
 * @author Brian, Raiden */
public final class NcrnaMitoSecondarySummary {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoSecondarySummary.class,false);
		try{final Parser parser=new Parser();String arm=null,verified=null;int expected=-1;
			for(String arg:pp.args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
				if(k.equals("arm")){arm=v;}else if(k.equals("expected")){expected=Parse.parseIntKMG(v);}else if(k.equals("verified")){verified=v;}
				else if(!parser.parse(arg,k,v)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			require(!parser.append && parser.maxReads<0 && verified!=null,"Use the complete saved member table and a producer-verified manifest, never append or a read limit");
			final ObjectMap<String,String> pins=EukPanelGffValidation.hashes(verified);require(parser.in1!=null && pins.get(parser.in1)!=null,"The saved member input must be explicitly pinned");
			summarize(parser.in1,parser.out1,arm,expected);
			pp.outstream.println("MITO_SECONDARY_SUMMARY_PASS arm="+arm+" copies="+expected+" strata=4 saved_identities=true no_alignments=true");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public static void summarize(String input,String output,String arm,int expected){
		require(input!=null && output!=null && arm!=null && !arm.isEmpty() && arm.indexOf('\t')<0 && arm.indexOf('\n')<0 && expected>0,"Require a saved input, fresh output, literal arm label and positive complete-copy denominator");
		RrnaResourceIO.fresh(output,"");
		final ObjectIntMap<String> loci=new ObjectIntMap<String>(String.class),ids=new ObjectIntMap<String>(String.class);
		final ArrayList<Unique> unique=new ArrayList<Unique>();final Summary[] sums={new Summary(),new Summary(),new Summary(),new Summary()};
		final ByteFile in=RrnaResourceIO.open(input);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,HEADER,input);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==9,"Saved secondary member schema changed");
				final String locus=p.parseString(0),uid=p.parseString(1);final int length=p.parseInt(2);final double identity=p.parseDouble(5);
				require(!locus.isEmpty() && !uid.isEmpty() && loci.get(locus)<0 && length>0,"Every original copy needs a unique nonempty locus and retained sequence ID/length");loci.put(locus,1);
				require(Double.isFinite(identity) && identity>=0 && identity<=1,"Saved member identity must be finite and in[0,1]; no missing or invented values");
				final boolean partial=flag(p.parseString(3)),pseudo=flag(p.parseString(4));int index=ids.get(uid);
				if(index<0){index=unique.size();ids.put(uid,index);unique.add(new Unique(length,identity));}
				else{final Unique u=unique.get(index);require(u.length==length && u.identity==identity,"All copies of a unique sequence shared one final alignment, so length and saved identity must agree: "+uid);}
				sums[0].add(index,identity);if(!partial && !pseudo){sums[1].add(index,identity);}if(partial){sums[2].add(index,identity);}if(pseudo){sums[3].add(index,identity);}
			}
		}finally{require(!in.close(),"Saved member table read failed");}
		require(loci.size()==expected && sums[0].copies==expected && sums[0].seen.size()==unique.size(),"Summary must conserve all original copies and exact unique IDs before writing output");
		final ByteStreamWriter out=RrnaResourceIO.writer(output);final ByteBuilder b=new ByteBuilder();
		try{out.println("arm\tstratum\tcopies\tunique_sequences\tcopy_mean_identity\tunique_mean_identity\tidentity_source");
			for(int i=0;i<sums.length;i++){final Summary s=sums[i];b.clear().append(arm).tab().append(STRATA[i]).tab().append(s.copies).tab().append(s.seen.size()).tab();
				if(s.copies==0){b.append("NA\tNA");}else{b.append(s.copySum/s.copies,9).tab().append(s.uniqueSum/s.seen.size(),9);}out.print(b.append("\tSAVED_NINE_DECIMAL_QUANTUM_IDENTITY").nl());}
		}finally{RrnaResourceIO.finish(out);}
	}
	private static boolean flag(String value){require(value.equals("true") || value.equals("false"),"Source partial/pseudo flags must remain explicit booleans");return value.equals("true");}
	private static void require(boolean value,String message){if(!value){throw new IllegalArgumentException(message);}}
	private static final class Unique{Unique(int n,double d){length=n;identity=d;}final int length;final double identity;}
	private static final class Summary{
		void add(int index,double identity){assert(index>=0):"Unique indices come from the validated global sequence roster";copies++;copySum+=identity;if(seen.get(index)<0){seen.put(index,1);uniqueSum+=identity;}}
		int copies;double copySum,uniqueSum;final IntHashMap seen=new IntHashMap();
	}
	public static final String HEADER="locus_id\tunique_id\tlength\tpartial\tpseudo\tquantum_identity\taligned_bases\tconsensus_start\tconsensus_stop";
	private static final String[] STRATA={"ALL","SOURCE_NONPARTIAL_NONPSEUDO","PARTIAL","PSEUDO"};
}
