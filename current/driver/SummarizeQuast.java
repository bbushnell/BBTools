package driver;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date May 8, 2015
 *
 */
public class SummarizeQuast {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Create a new SummarizeQuast instance
		SummarizeQuast sq=new SummarizeQuast(args);
		
		///And run it
		LinkedHashMap<String, LinkedHashMap<String, ArrayList<Double>>> map=sq.summarize();
		
		sq.print(map);
	}
	
	/**
	 * Constructs SummarizeQuast instance and parses command-line arguments.
	 * Processes input file paths, output destination, and analysis options.
	 * @param args Command line arguments including file paths and options
	 */
	public SummarizeQuast(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		ArrayList<String> names=new ArrayList<String>();
		Parser parser=new Parser();
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("required")){
				requiredString=b;
			}else if(a.equals("normalize")){
				normalize=Parse.parseBoolean(b);
			}else if(a.equals("box")){
				box=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(!arg.contains("=")){
				String[] x=(new File(arg).exists() ? new String[] {arg} : arg.split(","));
				for(String x2 : x){
					if(new File(x2).exists()){
						names.add(x2);
					}
				}
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		
		{//Process parser fields
			out=(parser.out1==null ? "stdout" : parser.out1);
			if(parser.in1!=null){
				String[] x=(new File(parser.in1).exists() ? new String[] {parser.in1} : parser.in1.split(","));
				for(String x2 : x){names.add(x2);}
			}
		}

		in=new ArrayList<String>();
		for(String s : names){
			Tools.getFileOrFiles(s, in, false, false, false, true);
		}
	}
	
	/**
	 * Main processing method that aggregates metrics from all QUAST files.
	 * Creates QuastSummary objects for each input file, optionally normalizes
	 * values, and consolidates metrics by name and assembly.
	 * @return Nested map: metric name -> assembly name -> list of values
	 */
	public LinkedHashMap<String, LinkedHashMap<String, ArrayList<Double>>> summarize(){
		
		ArrayList<QuastSummary> alqs=new ArrayList<QuastSummary>();
		for(String path : in){
			QuastSummary qs=new QuastSummary(path);
			if(normalize){qs.normalize();}
			alqs.add(qs);
		}

		LinkedHashMap<String, LinkedHashMap<String, ArrayList<Double>>> metricMap=new LinkedHashMap<String, LinkedHashMap<String, ArrayList<Double>>>();
		for(QuastSummary qs : alqs){
			for(String metricName : qs.metrics.keySet()){
				LinkedHashMap<String, ArrayList<Double>> asmMap=metricMap.get(metricName);
				if(asmMap==null){
					asmMap=new LinkedHashMap<String, ArrayList<Double>>();
					metricMap.put(metricName, asmMap);
				}
				ArrayList<Entry> ale=qs.metrics.get(metricName);
				assert(ale!=null);
//				assert(!ale.isEmpty()) : qs.path+"\n"+metricName+"\n";
				for(Entry e : ale){
					ArrayList<Double> ald=asmMap.get(e.assembly);
					if(ald==null){
						ald=new ArrayList<Double>();
						asmMap.put(e.assembly, ald);
					}
					ald.add(e.value);
				}
			}
		}
		
		return metricMap;
	}
	
	/**
	 * Outputs aggregated metrics to the specified destination.
	 * Formats results as tab-delimited text with metric names as headers.
	 * Optionally outputs box plot statistics (10th, 25th, 50th, 75th, 90th percentiles)
	 * instead of individual values when box mode is enabled.
	 *
	 * @param metricMap Nested map containing aggregated metric values
	 */
	public void print(LinkedHashMap<String, LinkedHashMap<String, ArrayList<Double>>> metricMap){
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		tsw.start();
		for(String metricName : metricMap.keySet()){
			LinkedHashMap<String, ArrayList<Double>> asmMap=metricMap.get(metricName);
			if(asmMap!=null && !asmMap.isEmpty()){
				tsw.println("\n"+metricName);
				assert(!asmMap.isEmpty());
				for(String asm : asmMap.keySet()){
					ArrayList<Double> ald=asmMap.get(asm);
					assert(ald!=null);
					assert(!ald.isEmpty());
					if(ald!=null && !ald.isEmpty()){
						tsw.print(asm);
						if(box){
							double[] array=new double[ald.size()];
							for(int i=0; i<ald.size(); i++){array[i]=ald.get(i);}
							Arrays.sort(array);
							final int len=array.length-1;
							tsw.print("\t"+array[(int)Math.round(0.1*len)]);
							tsw.print("\t"+array[(int)Math.round(0.25*len)]);
							tsw.print("\t"+array[(int)Math.round(0.5*len)]);
							tsw.print("\t"+array[(int)Math.round(0.75*len)]);
							tsw.print("\t"+array[(int)Math.round(0.9*len)]);
						}else{
							for(Double d : ald){
								tsw.print("\t"+d);
							}
						}
						tsw.println();
					}
				}
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Represents metrics extracted from a single QUAST output file.
	 * Parses tabular QUAST format and stores metrics organized by row name.
	 * Supports value normalization by dividing each metric by its mean.
	 */
	private class QuastSummary{
		
		/** Constructs QuastSummary by parsing the specified QUAST file.
		 * @param path_ Path to the QUAST output file */
		QuastSummary(String path_){
			path=path_;
			metrics=process(path);
		}
		
		/**
		 * Parses a QUAST tabular file and extracts metrics.
		 * Reads header row to identify assembly columns, then processes data rows
		 * to create Entry objects for each numeric metric value.
		 * Filters assemblies based on requiredString if specified.
		 *
		 * @param fname Path to the QUAST file to parse
		 * @return Map of metric name to list of Entry objects
		 */
		LinkedHashMap<String, ArrayList<Entry>> process(String fname){
			LinkedHashMap<String, ArrayList<Entry>> map=new LinkedHashMap<String, ArrayList<Entry>>();
			TextFile tf=new TextFile(fname);
			String[] header=null;
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				String[] split=line.split("\t");
				if(header==null){
					header=split;
				}else{
					final String row=split[0];
					ArrayList<Entry> list=new ArrayList<Entry>(split.length-1);
					map.put(row, list);
					for(int i=1; i<split.length; i++){
						final String col=header[i];
						if(requiredString==null || col.contains(requiredString)){
							try {
								Entry e=new Entry(col, split[i]);
								if(!Double.isNaN(e.value) && !Double.isInfinite(e.value)){
									list.add(e);
								}
							} catch (NumberFormatException ex) {
								//Do nothing
							}
						}
					}
					
//					assert(false) : row+", "+list+", "+split.length+"\n"+Arrays.toString(split);
					
				}
			}
			return map;
		}
		
		/** Normalizes all metric values by dividing each by the mean of its metric group.
		 * Applied to make metrics with different scales comparable. */
		void normalize(){
			for(ArrayList<Entry> list : metrics.values()){
				normalize(list);
			}
		}
		
		/**
		 * Normalizes a single metric's values by dividing each by the group mean.
		 * Handles empty lists and zero averages to avoid division errors.
		 * @param list List of Entry objects to normalize
		 */
		private void normalize(ArrayList<Entry> list){
			if(list.isEmpty()){return;}
			if(list==null || list.isEmpty()){return;}
			double sum=0;
			for(Entry e : list){
				sum+=e.value;
			}
			double avg=sum/list.size();
			double mult=(avg==0 ? 1 : 1/avg);
			for(Entry e : list){
				e.value*=mult;
			}
		}
		
		/** Map of metric names to lists of Entry objects containing assembly values */
		final LinkedHashMap<String, ArrayList<Entry>> metrics;
		/** Path to the source QUAST file */
		final String path;
		
	}
	
	/** Represents a single metric value for a specific assembly.
	 * Pairs an assembly name with its corresponding numeric metric value. */
	private class Entry{
		
		/**
		 * Constructs Entry by parsing string value as double.
		 * @param assembly_ Assembly name identifier
		 * @param value_ String representation of metric value
		 * @throws NumberFormatException if value cannot be parsed as double
		 */
		Entry(String assembly_, String value_) throws NumberFormatException {
			this(assembly_, Double.parseDouble(value_));
		}
		
		/**
		 * Constructs Entry with pre-parsed numeric value.
		 * @param assembly_ Assembly name identifier
		 * @param value_ Numeric metric value
		 */
		Entry(String assembly_, double value_){
			assembly=assembly_;
			value=value_;
		}
		
		/** Assembly name identifier */
		String assembly;
		/** Numeric metric value for this assembly */
		double value;
		
	}
	
	/** List of input QUAST file paths to process */
	final ArrayList<String> in;
	/** Output destination for aggregated results */
	final String out;
	
	/**
	 * Optional filter string - only assemblies containing this string are included
	 */
	String requiredString=null;
	/** Whether to normalize metric values by dividing by group means */
	boolean normalize=true;
	/** Whether to output box plot percentiles instead of individual values */
	boolean box=true;
	
}
