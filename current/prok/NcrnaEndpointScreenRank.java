package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Ranks already measured native endpoint-count rows; it never learns tables or
 * networks. Uses the text template's Parser, byte reader and owned output pattern.
 * Retained geometry records own primitive counts. Exact fractions avoid rounding
 * ties; equal placements remain tied rather than selecting an invented offset.
 * VALIDATION decides ordering; TRAIN_LOO is displayed separately. The result is a
 * count-screening shortlist, not significance, neural accuracy or caller adoption.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaEndpointScreenRank{

	NcrnaEndpointScreenRank(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("mode")){check("rank".equals(b), "Saved-row ranking needs mode=rank");}
			else{check(parser.parse(arg, a, b), "Unknown endpoint-ranking argument: "+arg);}
		}
		in=parser.in1;
		out=parser.out1;
		check(in!=null && out!=null && parser.in2==null && parser.out2==null && !parser.overwrite
			&& !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Supply mode=rank in=native_curves.tsv out=fresh_prefix without filtering or overwrite");
		check(Tools.testInputFiles(false, true, in)
			&& Tools.testOutputFiles(false, false, false, out+".ranked.tsv", out+".choices.tsv", out+".byk.tsv"), "Unreadable curves or existing ranking outputs");
		Shared.setThreads(1);
	}

	void process(){
		final Geometry[][] geometries=read();
		try(Output ranked=new Output(out+".ranked.tsv"); Output choices=new Output(out+".choices.tsv");
				Output byK=new Output(out+".byk.tsv")){
			ranked.writer.println(HEADER);
			choices.writer.println(HEADER);
			byK.writer.println(HEADER);
			for(int end=0; end<2; end++){
				final ArrayList<Geometry> present=new ArrayList<Geometry>();
				for(final Geometry geometry:geometries[end]){if(geometry!=null){present.add(geometry);}}
				final Geometry[] order=present.toArray(new Geometry[present.size()]);
				Arrays.sort(order, (a, b)->{
					final int result=compare(a, b);
					return result!=0 ? result : Integer.compare(a.outside, b.outside);
				});
				final Geometry[] bestByK=new Geometry[16];
				int rank=0;
				boolean available=false;
				Geometry previous=null;
				for(int index=0; index<order.length; index++){
					final Geometry geometry=order[index];
					final boolean valid=geometry.n[1]>0, top=valid && compare(order[0], geometry)==0;
					if(bestByK[geometry.k]==null){bestByK[geometry.k]=geometry;}
					if(valid && (previous==null || compare(previous, geometry)!=0)){rank=index+1;}
					write(ranked, geometry, valid ? rank : 0, top ? "TOP_OR_TIED_COUNT_SCREEN" : valid ? "RANKED_COUNT_SCREEN" : "UNAVAILABLE_NO_VALIDATION");
					if(top){
						write(choices, geometry, rank, "TOP_OR_TIED_COUNT_SCREEN");
						available=true;
					}
					previous=geometry;
				}
				for(int k=1; k<bestByK.length; k++){
					for(final Geometry geometry:order){
						if(geometry.k!=k || compare(bestByK[k], geometry)!=0){continue;}
						final boolean valid=geometry.n[1]>0;
						write(byK, geometry, valid ? 1 : 0, valid ? "TOP_OR_TIED_AT_K_COUNT_SCREEN" : "UNAVAILABLE_NO_VALIDATION");
					}
				}
				if(!available){
					check(order.length>0, "Both native endpoints must remain in ranking accounting");
					choices.writer.print(buffer.clear().append(end==0 ? "START" : "STOP")
						.append("\tNA\tNA\tNA\tNA\tNA\t0\t0\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tUNAVAILABLE_NO_VALIDATION\tCOUNTS_NOT_SIGNIFICANCE_NEURAL_ACCURACY_OR_CALLER_SELECTION\tNA\tNA\n"));
				}
			}
		}
	}

	/** Native populations are int-indexed, so exact count products fit signed long. */
	private Geometry[][] read(){
		final Geometry[][] result=new Geometry[2][136];
		final LineParser1 p=new LineParser1('\t');
		int measuredRadius=-1;
		try(Input input=new Input(in)){
			RrnaResourceIO.header(input.reader, CURVES, in);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==13, "Native endpoint-curve schema changed");
				final int k=number(p, 0), inside=number(p, 1), outside=number(p, 2), radius=number(p, 5);
				final int end=p.termEquals("START", 3) ? 0 : p.termEquals("STOP", 3) ? 1 : -1;
				final int phase=p.termEquals("TRAIN_LOO", 4) ? 0 : p.termEquals("VALIDATION", 4) ? 1 : -1;
				check(k>=1 && k<=15 && inside+outside==k && end>=0 && phase>=0 && radius>0,
					"Ranking needs native semantic geometry and separate TRAIN_LOO/VALIDATION rows");
				if(measuredRadius<0){measuredRadius=radius;}
				check(radius==measuredRadius, "A saved native screen has one fixed candidate radius");
				final int key=k*(k+1)/2+outside;
				Geometry geometry=result[end][key];
				if(geometry==null){
					geometry=new Geometry(end, k, inside, outside, radius);
					result[end][key]=geometry;
				}
				// Nonranking rows still declare a geometry that must have all target rows.
				final int ratio=p.termEquals("2", 6) || p.termEquals("2.0", 6) ? 0
					: p.termEquals("3", 6) || p.termEquals("3.0", 6) ? 1 : -1;
				if(ratio<0){continue;}
				final int bit=1<<(phase*2+ratio);
				check((geometry.seen&bit)==0, "Duplicate target ratio for one geometry");
				final long n=number(p, 7), correct=number(p, 8), other=number(p, 9), abstain=number(p, 10);
				check(correct+other+abstain==n, "Native count outcomes must conserve their int-indexed eligible population");
				if((geometry.seen&(3<<(phase*2)))!=0){check(geometry.n[phase]==n, "Ratios two and three must use the same native cohort");}
				geometry.n[phase]=n;
				if(ratio==0){geometry.two[phase]=correct;}
				else{geometry.three[phase]=correct;}
				geometry.seen|=bit;
			}
		}
		for(int end=0; end<2; end++){
			int count=0;
			for(final Geometry geometry:result[end]){
				if(geometry==null){continue;}
				check(geometry.seen==15 && geometry.three[0]<=geometry.two[0] && geometry.three[1]<=geometry.two[1],
					"Each geometry needs both ratios/partitions and monotone decisive-true counts");
				count++;
			}
			check(count>0, "Both endpoints need native ratios two and three before ranking");
		}
		return result;
	}

	/** Selection order excludes outside placement; its display ordering never breaks a tie. */
	private static int compare(final Geometry a, final Geometry b){
		if(a.n[1]==0 || b.n[1]==0){return a.n[1]==b.n[1] ? Integer.compare(a.k, b.k) : a.n[1]==0 ? 1 : -1;}
		int result=Long.compare(b.two[1]*a.n[1], a.two[1]*b.n[1]);
		if(result==0){result=Long.compare(b.three[1]*a.n[1], a.three[1]*b.n[1]);}
		return result==0 ? Integer.compare(a.k, b.k) : result;
	}

	private void write(final Output output, final Geometry geometry, final int rank, final String status){
		assert(geometry.seen==15) : "Only complete native ratio/partition records may enter the saved ranking";
		buffer.clear().append(geometry.end==0 ? "START" : "STOP").tab();
		if(rank==0){buffer.append("NA");}else{buffer.append(rank);}
		buffer.tab().append(geometry.k).tab().append(geometry.inside).tab().append(geometry.outside).tab().append(geometry.radius);
		for(int slot=0; slot<2; slot++){
			final int phase=1-slot;
			buffer.tab().append(geometry.n[phase]).tab().append(geometry.two[phase]);
			fraction(geometry.two[phase], geometry.n[phase]);
			buffer.tab().append(geometry.three[phase]);
			fraction(geometry.three[phase], geometry.n[phase]);
		}
		buffer.tab().append(status).append("\tCOUNTS_NOT_SIGNIFICANCE_NEURAL_ACCURACY_OR_CALLER_SELECTION");
		gap(geometry.two, geometry.n);
		gap(geometry.three, geometry.n);
		output.writer.print(buffer.nl());
	}
	/** Positive means a larger training fraction; unavailable cohorts stay explicit. */
	private void gap(final long[] numerator, final long[] denominator){
		assert(numerator.length==2 && denominator.length==2) : "Gap compares TRAIN_LOO and VALIDATION on their own eligible cohorts";
		buffer.tab();
		if(denominator[0]==0 || denominator[1]==0){buffer.append("NA");}
		else{buffer.append(numerator[0]/(double)denominator[0]-numerator[1]/(double)denominator[1], 8);}
	}
	private void fraction(final long numerator, final long denominator){
		assert(numerator>=0 && numerator<=denominator) : "Published fractions must use the validated native eligible cohort";
		buffer.tab();
		if(denominator==0){buffer.append("NA");}else{buffer.append(numerator/(double)denominator, 8);}
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static int number(final LineParser1 parser, final int term){
		parser.setBounds(term);
		return NcrnaNumbers.unsignedInt(parser.line(), parser.a(), parser.b());
	}

	private static final class Geometry{
		Geometry(final int end_, final int k_, final int inside_, final int outside_, final int radius_){
			end=end_;
			k=k_;
			inside=inside_;
			outside=outside_;
			radius=radius_;
		}
		final int end, k, inside, outside, radius;
		final long[] n=new long[2], two=new long[2], three=new long[2];
		int seen;
	}
	private final String in, out;
	private final ByteBuilder buffer=new ByteBuilder();
	private static final String CURVES="k\tinside\toutside\tboundary\tpartition\tradius\tratio\teligible_loci\tdecisive_true\tdecisive_other\tabstain\ttrue_fraction\tother_fraction";
	private static final String HEADER="boundary\trank\tk\tinside\toutside\tradius\tvalidation_n\tvalidation_true2\tvalidation_fraction2\tvalidation_true3\tvalidation_fraction3\ttrain_n\ttrain_true2\ttrain_fraction2\ttrain_true3\ttrain_fraction3\tstatus\tqualification\ttrain_minus_validation2\ttrain_minus_validation3";
}
