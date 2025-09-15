package fun;

import java.util.HashMap;

public class Merced {

	public static void main(String[] args) {
		
		//Create a new Map, mapping Integers to Strings
		HashMap<String, Integer> map=new HashMap<String, Integer>();
		//Associate value 21 with key "Age"
//		map.put("Age", 21);
//		map.put("Weight", 130);
//		//This won't compile since 17.0 is a Double, not an int
//		map.put("Age", 17.0);
//		map.put("Age", 23);
//		//This won't compile since "Tall" is a String, not int
//		map.put("Height", "Tall");
//		//What will this return?
//		Integer age=map.get("Age");
//		//This will succeed, but return null ("None" in Python)
//		Integer ssn=map.get("SSN#");
//		//Won't compile because the return is an Integer
//		String weight=map.get("Weight");
	}
	
	class Person {
		
		Person(String name_){
			name=name_;
		}
		
		float attractiveness() {
			return looks*salary/age;
		}
		
		String name;
		int age;
		int salary;
		float looks;
		float brains;//TODO: Irrelevant, remove
		
	}
	
	HashMap<String, Person> nameTable=new HashMap<String, Person>();
	
}
