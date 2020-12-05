// test runner

use crate::algorithm::relate::IntersectionMatrix;
use geo_types::Geometry;

#[derive(Debug)]
struct RelateTestFailure {
    actual_result: IntersectionMatrix,
    input: RelateTestCase,
}

#[derive(Debug)]
struct TestRunner {
    dir_path: String,
    failures: Vec<RelateTestFailure>,
    successes: Vec<RelateTestCase>,
}

#[derive(Debug)]
struct RelateTestCase {
    expected_result: IntersectionMatrix,
    test_file_name: String,
    description: String,
    geometry_a: Geometry<f64>,
    geometry_b: Geometry<f64>,
}

/// <run>
/// <precisionModel scale="1.0" offsetx="0.0" offsety="0.0"/>
///
/// <case>
/// <desc>AA disjoint</desc>
/// <a>
/// POLYGON(
/// (0 0, 80 0, 80 80, 0 80, 0 0))
/// </a>
/// <b>
/// POLYGON(
/// (100 200, 100 140, 180 140, 180 200, 100 200))
/// </b>
/// <test><op name="relate" arg3="FF2FF1212" arg1="A" arg2="B"> true </op>
/// </test>
/// <test>  <op name="intersects" arg1="A" arg2="B">   false   </op></test>
/// <test>  <op name="contains" arg1="A" arg2="B">   false   </op></test>
/// </case>
/// </run>
mod test_input {
    use geo_types::Geometry;
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    pub(crate) struct Run {
        // TODO?
        // precision_model: PrecisionModel
        #[serde(rename = "case")]
        pub cases: Vec<Case>,
    }

    #[derive(Debug, Deserialize)]
    pub(crate) struct Case {
        pub(crate) desc: String,
        pub(crate) a: String, // WKT
        pub(crate) b: String, // WKT
        #[serde(rename = "test")]
        pub(crate) tests: Vec<Test>,
    }

    #[derive(Debug, Deserialize)]
    pub(crate) struct Test {
        pub(crate) op: Operation,
    }

    impl Test {
        pub(crate) fn arg1(&self) -> Option<&str> {
            todo!()
        }
        pub(crate) fn arg2(&self) -> Option<&str> {
            todo!()
        }
        pub(crate) fn arg3(&self) -> Option<&str> {
            todo!()
        }
    }

    #[derive(Debug, Deserialize)]
    pub(crate) struct Operation {
        pub(crate) name: String,
        // seems to always be "A" or "B", indicating the geometry in the containing Case
        pub(crate) arg1: String,
        pub(crate) arg2: String,

        // If set, used to pass in the expected DENIM matrix
        pub(crate) arg3: Option<String>,

        // TODO: not yet using this. Sometimes a bool, sometimes wkt, depends on the operation
        #[serde(rename = "$value")]
        pub(crate) expected: String,
    }
}

impl TestRunner {
    fn new(dir_path: String) -> Self {
        TestRunner {
            dir_path,
            successes: Vec::new(),
            failures: Vec::new(),
        }
    }

    #[track_caller]
    fn run_all(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let cases = self.parse_cases()?;
        println!("cases.len(): {}", cases.len());
        for case in cases {
            use crate::algorithm::relate::relate_computer::RelateComputer;

            let mut relate_computer = RelateComputer::new(&case.geometry_a, &case.geometry_b);
            let intersection_matrix = relate_computer.compute_intersection_matrix();
            if intersection_matrix == case.expected_result {
                println!("succeeded: {:?}", case);
                self.successes.push(case);
            } else {
                let failure = RelateTestFailure {
                    input: case,
                    actual_result: intersection_matrix,
                };
                println!("failed: {:?}", failure);
                self.failures.push(failure);
            }
        }
        Ok(())
    }

    #[track_caller]
    fn parse_cases(&self) -> Result<Vec<RelateTestCase>, Box<dyn std::error::Error>> {
        use std::convert::TryFrom;
        use std::fs::File;
        use std::io::BufReader;

        let mut cases = vec![];
        for entry in std::fs::read_dir(&self.dir_path)? {
            let entry = entry?;
            if !entry
                .file_name()
                .to_str()
                .map(|name| name.contains("Relate"))
                .unwrap_or(false)
            {
                // ignoring non-Relate file
                continue;
            }

            let file = File::open(entry.path())?;
            let file_reader = BufReader::new(file);
            let run: test_input::Run = serde_xml_rs::from_reader(file_reader).map_err(|err| {
                format!("invalid test input: {:?}. error: {:?}", entry.path(), err)
            })?;
            for case in run.cases {
                // re: `unwrap` see https://github.com/georust/wkt/issues/49
                // we could map to our own error or something, but this is just test code anyway
                let wkt_a = match wkt::Wkt::from_str(&case.a) {
                    Ok(wkt) => wkt,
                    Err(e) => {
                        println!("error: {:?}, skipping invalid WKT: {}", e, &case.a);
                        continue;
                    }
                };
                let geometry_a = Geometry::try_from(wkt_a).expect("conversion error");

                let wkt_b = match wkt::Wkt::from_str(&case.b) {
                    Ok(wkt) => wkt,
                    Err(e) => {
                        println!("error: {:?}, skipping invalid WKT: {}", e, &case.b);
                        continue;
                    }
                };
                let geometry_b = Geometry::try_from(wkt_b).expect("conversion error");

                for test in case.tests {
                    match test.op.name.as_str() {
                        "relate" => {
                            let expected_result = match (
                                test.op.arg1.as_str(),
                                test.op.arg2.as_str(),
                                test.op.arg3.as_ref(),
                            ) {
                                ("A", "B", Some(im_str)) => IntersectionMatrix::from_str(im_str),
                                _ => panic!("unexpected args for relate operation"),
                            };
                            let relate_test_case = RelateTestCase {
                                test_file_name: entry
                                    .path()
                                    .file_name()
                                    .unwrap()
                                    .to_string_lossy()
                                    .to_string(),
                                description: case.desc.clone(),
                                geometry_a: geometry_a.clone(),
                                geometry_b: geometry_b.clone(),
                                expected_result,
                            };
                            cases.push(relate_test_case);
                        }
                        _ => {} // ignoring non-relate test cases for now
                    }
                }
            }
        }
        Ok(cases)
    }
}

#[test]
fn test_general_cases() {
    let mut runner =
        TestRunner::new("/Users/mkirk/src/georust/geo/geo/resources/testxml/general".to_string());
    runner.run_all().expect("error while running tests");
    assert!(
        !runner.successes.is_empty(),
        "successes: {:?}",
        &runner.successes
    );

    let total = runner.successes.len() + runner.failures.len();
    assert!(
        runner.failures.is_empty(),
        "test runner failed {} of {} cases. Failures: {:?}",
        runner.failures.len(),
        total,
        &runner.failures
    );
}
