package fileIO;

import java.util.ArrayList;

/**
 * TextFile subclass that collects remaining logical lines into a String array.
 * Uses the inherited blank-line predicate; its additional filter placeholder is
 * inactive. Instances retain TextFile's mutable reader state, and callers own closure.
 * @author Brian Bushnell
 */
public class GenericTextFile extends TextFile{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens the input immediately with subprocess permission disabled.
	 * ReadWrite still selects the actual backend, as for TextFile.
	 * @param name Input name interpreted by ReadWrite */
	public GenericTextFile(String name){
		super(name, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Collects remaining logical lines in order through virtual nextLine() calls.
	 * Stops at the first null result without resetting or closing the reader.
	 * Subclasses overriding nextLine() can change which lines are collected.
	 * @return New array of retained lines, empty when none remain
	 */
	public String[] toLines(){
		String s=null;
		ArrayList<String> list=new ArrayList<String>(4096);
		for(s=nextLine(); s!=null; s=nextLine()){
			list.add(s);
		}
		return list.toArray(new String[list.size()]);
	}

	/** Returns the next TextFile logical line; the extra filter loop is inactive.
	 * Retained lines are not trimmed, and EOF does not close the reader.
	 * @return Next retained line, or null at EOF or when the reader is closed */
	@Override
	public String nextLine(){
		String line=readLine();
		while(line!=null && false){
			line=readLine();
		}
		return line;
	}
}
