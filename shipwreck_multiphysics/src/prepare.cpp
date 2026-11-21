#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;
using Point = std::array<double, 3>;
using Triangle = std::array<Point, 3>;

struct Flow {
    double time, speed, concentration;
};

struct Inflow {
    std::vector<Flow> rows;
    std::string csv;
};

double number(const std::string& text) {
    std::size_t end = 0;
    const double value = std::stod(text, &end);
    if (end != text.size() || !std::isfinite(value)) throw std::runtime_error("Invalid number: " + text);
    return value;
}

Inflow readFlow(const fs::path& path, double endTime) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open flow CSV: " + path.string());

    Inflow inflow;
    std::string line;
    auto nextRow = [&] {
        if (!std::getline(input, line)) return false;
        // Preserve the exact input, including CRLF and a missing final newline.
        inflow.csv.append(line);
        if (!input.eof()) inflow.csv.push_back('\n');
        if (!line.empty() && line.back() == '\r') line.pop_back();
        return true;
    };
    if (!nextRow()) throw std::runtime_error("Empty flow CSV");
    if (line != "time_s,speed_m_s,Cs_volume_fraction") throw std::runtime_error("Invalid flow CSV header");

    auto& rows = inflow.rows;
    while (nextRow()) {
        const auto first = line.find(',');
        const auto second = first == std::string::npos ? first : line.find(',', first+1);
        const auto row = std::to_string(rows.size()+2);
        if (second == std::string::npos || line.find(',', second+1) != std::string::npos) {
            throw std::runtime_error("Expected three comma-separated values at flow row " + row);
        }
        auto readValue = [&](std::size_t position, std::size_t length = std::string::npos) {
            try {
                return number(line.substr(position, length));
            } catch (const std::exception& error) {
                throw std::runtime_error("Invalid flow row " + row + ": " + error.what());
            }
        };
        const Flow f{readValue(0, first), readValue(first+1, second-first-1), readValue(second+1)};
        if (f.speed <= 0 || f.concentration < 0 || f.concentration > .1 ||
            (rows.empty() ? f.time != 0 : f.time <= rows.back().time)) {
            throw std::runtime_error("Invalid flow row " + row);
        }
        rows.push_back(f);
    }
    if (input.bad() || !input.eof()) throw std::runtime_error("Could not read flow CSV: " + path.string());
    if (rows.size() < 2 || rows.back().time < endTime) throw std::runtime_error("Flow must cover 0 through end time");
    return inflow;
}

Point subtract(const Point& a, const Point& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

Point cross(const Point& a, const Point& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}

double dot(const Point& a, const Point& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

std::uint32_t little32(const unsigned char* bytes) {
    return std::uint32_t(bytes[0]) | std::uint32_t(bytes[1]) << 8 |
           std::uint32_t(bytes[2]) << 16 | std::uint32_t(bytes[3]) << 24;
}

std::vector<Triangle> readHull(const fs::path& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    const auto fileSize = static_cast<std::streamoff>(input.tellg());
    input.seekg(0);
    std::array<unsigned char, 84> header{};
    if (!input.read(reinterpret_cast<char*>(header.data()), header.size())) throw std::runtime_error("Cannot read binary STL");
    const auto count = little32(header.data() + 80);
    if (count == 0 || count > 2000000 || fileSize != static_cast<std::streamoff>(84ULL + 50ULL * count)) {
        throw std::runtime_error("Expected a complete binary STL with 1 to 2,000,000 triangles");
    }
    std::vector<Triangle> faces;
    faces.reserve(count);
    for (std::uint32_t i = 0; i < count; ++i) {
        std::array<unsigned char, 50> record{};
        if (!input.read(reinterpret_cast<char*>(record.data()), record.size())) throw std::runtime_error("Truncated STL");
        Triangle face{};
        for (std::size_t j = 0; j < 9; ++j) {
            const auto bits = little32(record.data() + 12 + 4*j);
            float value;
            static_assert(sizeof(value) == sizeof(bits) && std::numeric_limits<float>::is_iec559);
            std::memcpy(&value, &bits, sizeof(value));
            if (!std::isfinite(value)) throw std::runtime_error("Nonfinite STL vertex");
            face[j/3][j%3] = value;
        }
        faces.push_back(face);
    }
    return faces;
}

std::array<Point,2> bounds(const std::vector<Triangle>& faces) {
    std::array<Point,2> box{faces.front().front(), faces.front().front()};
    for (const auto& face : faces) for (const auto& p : face) for (std::size_t j = 0; j < 3; ++j) {
        box[0][j] = std::min(box[0][j], p[j]);
        box[1][j] = std::max(box[1][j], p[j]);
    }
    return box;
}

void verifyClosed(std::vector<Triangle>& faces) {
    using Edge = std::pair<Point, Point>;
    struct Edges { unsigned count = 0; int orientation = 0; };
    std::map<Edge, Edges> edges;
    double sixVolume = 0;
    const auto origin = faces.front().front();
    for (const auto& f : faces) {
        const auto normal = cross(subtract(f[1], f[0]), subtract(f[2], f[0]));
        const double normSquared = dot(normal, normal);
        if (!(normSquared > 0) || !std::isfinite(normSquared)) {
            throw std::runtime_error("STL triangle must have a finite, nonzero normal");
        }
        sixVolume += dot(subtract(f[0], origin), cross(subtract(f[1], origin), subtract(f[2], origin)));
        for (int i = 0; i < 3; ++i) {
            const auto& a = f[i];
            const auto& b = f[(i+1)%3];
            auto& edge = edges[std::minmax(a, b)];
            ++edge.count;
            edge.orientation += a < b ? 1 : -1;
        }
    }
    for (const auto& [edge, use] : edges) {
        if (use.count != 2 || use.orientation != 0) throw std::runtime_error("STL edges must form a closed, consistently wound surface");
    }
    if (!std::isfinite(sixVolume) || sixVolume == 0) throw std::runtime_error("STL has no finite nonzero signed volume");
    if (sixVolume < 0) for (auto& face : faces) std::swap(face[1], face[2]);
}

std::string tuple(const Point& p) {
    if (std::any_of(p.begin(), p.end(), [](double value){ return !std::isfinite(value); })) {
        throw std::runtime_error("Cannot write a nonfinite case vector");
    }
    std::ostringstream out;
    out << std::setprecision(16) << '(' << p[0] << ' ' << p[1] << ' ' << p[2] << ')';
    return out.str();
}

void write(const fs::path& path, const std::string& value) {
    fs::create_directories(path.parent_path());
    std::ofstream out(path, std::ios::binary);
    out << value;
    out.close();
    if (!out) throw std::runtime_error("Could not write " + path.string());
}

std::string dict(const std::string& object, const std::string& cls = "dictionary") {
    return "FoamFile { version 2.0; format ascii; class " + cls + "; object " + object + "; }\n";
}

std::string real(double value) {
    if (!std::isfinite(value)) throw std::runtime_error("Cannot write a nonfinite case value");
    std::ostringstream out;
    out << std::setprecision(16) << value;
    return out.str();
}

void writeHull(const fs::path& path, const std::vector<Triangle>& faces) {
    fs::create_directories(path.parent_path());
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Could not open " + path.string());
    out << std::scientific << std::setprecision(16) << "solid shipwreck\n";
    for (const auto& f : faces) {
        const auto n = cross(subtract(f[1], f[0]), subtract(f[2], f[0]));
        const double norm = std::sqrt(dot(n, n));
        out << "facet normal " << n[0]/norm << ' ' << n[1]/norm << ' ' << n[2]/norm << "\nouter loop\n";
        for (const auto& p : f) out << "vertex " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
        out << "endloop\nendfacet\n";
    }
    out << "endsolid shipwreck\n";
    out.close();
    if (!out) throw std::runtime_error("Could not write " + path.string());
}

std::string mesh(double span, double depth, int refine) {
    // Fixed extents for a size family, independent of yaw and placement.
    const std::array<double,4> x{-3*span, -.85*span, .85*span, 5*span};
    const std::array<double,4> y{-2*span, -.85*span, .85*span, 2*span};
    const std::array<int,3> nx{14,32,20}, ny{10,32,10};
    auto id = [](int i, int j, int k) { return i + 4*j + 16*k; };
    auto face = [](int a, int b, int c, int d) {
        return '('+std::to_string(a)+' '+std::to_string(b)+' '+std::to_string(c)+' '+std::to_string(d)+')';
    };
    std::ostringstream out;
    out << dict("blockMeshDict") << "scale 1;\nvertices (\n";
    for (double z : {0., depth}) for (double yy : y) for (double xx : x) out << tuple({xx,yy,z}) << '\n';
    out << ");\nblocks (\n";
    for (int j=0; j<3; ++j) for (int i=0; i<3; ++i) {
        out << "hex (" << id(i,j,0) << ' ' << id(i+1,j,0) << ' ' << id(i+1,j+1,0) << ' ' << id(i,j+1,0) << ' '
            << id(i,j,1) << ' ' << id(i+1,j,1) << ' ' << id(i+1,j+1,1) << ' ' << id(i,j+1,1) << ") ("
            << nx[i]*refine << ' ' << ny[j]*refine << ' ' << 28*refine << ") simpleGrading ("
            << (i==0?.3:i==2?4.:1.) << ' ' << (j==0?.45:j==2?2.2:1.) << " 8)\n";
    }
    out << ");\nedges ();\nboundary (\n";
    auto patch = [&](const char* name, const char* type, auto faces) {
        out << name << " { type " << type << "; faces (\n";
        faces();
        out << "); }\n";
    };
    patch("inlet","patch",[&]{for(int j=0;j<3;++j)out<<face(id(0,j,0),id(0,j,1),id(0,j+1,1),id(0,j+1,0))<<'\n';});
    patch("outlet","patch",[&]{for(int j=0;j<3;++j)out<<face(id(3,j,0),id(3,j+1,0),id(3,j+1,1),id(3,j,1))<<'\n';});
    patch("sides","symmetry",[&]{for(int i=0;i<3;++i)out<<face(id(i,0,0),id(i+1,0,0),id(i+1,0,1),id(i,0,1))<<'\n'<<face(id(i,3,0),id(i,3,1),id(i+1,3,1),id(i+1,3,0))<<'\n';});
    patch("bed","wall",[&]{for(int j=0;j<3;++j)for(int i=0;i<3;++i)out<<face(id(i,j,0),id(i,j+1,0),id(i+1,j+1,0),id(i+1,j,0))<<'\n';});
    patch("top","symmetryPlane",[&]{for(int j=0;j<3;++j)for(int i=0;i<3;++i)out<<face(id(i,j,1),id(i+1,j,1),id(i+1,j+1,1),id(i,j+1,1))<<'\n';});
    out << ");\nmergePatchPairs ();\n";
    return out.str();
}

template<class F>
std::string table(const std::vector<Flow>& flow, F value) {
    std::ostringstream out;
    out << "type uniformFixedValue; uniformValue { type table; outOfBounds clamp; interpolationScheme linear; values (\n";
    for (const auto& row : flow) out << '(' << real(row.time) << ' ' << value(row) << ")\n";
    out << "); } value uniform " << value(flow.front()) << ';';
    return out.str();
}

void volumeField(const fs::path& casePath, const std::string& name, const std::string& dims,
                 const std::string& initial, const std::string& inlet, const std::string& outlet,
                 const std::string& bed, bool vector = false) {
    write(casePath/"0"/name, dict(name, vector?"volVectorField":"volScalarField") +
        "dimensions " + dims + ";\ninternalField uniform " + initial + ";\nboundaryField\n{\n"
        "inlet {" + inlet + "}\noutlet {" + outlet + "}\nbed {" + bed + "}\n"
        "sides {type symmetry;}\ntop {type symmetryPlane;}\n}\n");
}

class StagedCase {
    fs::path path_;
    bool published_ = false;

public:
    explicit StagedCase(const fs::path& output) : path_(output) {
        fs::create_directories(output.parent_path());
        path_ += ".preparing-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
        if (!fs::create_directory(path_)) throw std::runtime_error("Cannot create staged case");
    }

    StagedCase(const StagedCase&) = delete;
    StagedCase& operator=(const StagedCase&) = delete;

    ~StagedCase() {
        if (!published_) {
            std::error_code ignored;
            fs::remove_all(path_, ignored);
        }
    }

    const fs::path& path() const { return path_; }

    void publish(const fs::path& output) {
        if (fs::exists(output) || fs::is_symlink(output)) throw std::runtime_error("Case already exists");
        fs::rename(path_, output);
        published_ = true;
    }
};

void prepare(const std::map<std::string,std::string>& arguments) {
    auto get = [&](const char* key, const std::string& fallback) {
        const auto i = arguments.find(key);
        return i == arguments.end() ? fallback : i->second;
    };
    auto value = [&](const char* key, const char* fallback) { return number(get(key, fallback)); };
    const fs::path output = fs::absolute(arguments.at("--case")).lexically_normal();
    const fs::path templatePath = arguments.at("--template");
    if (fs::exists(output) || fs::is_symlink(output)) throw std::runtime_error("Case already exists");
    const double length = number(arguments.at("--length"));
    const double yaw = value("--yaw", "0"), burial = value("--burial", "0");
    const double cx = value("--x", "0"), cy = value("--y", "0");
    const double depth = number(arguments.at("--depth")), endTime = number(arguments.at("--end"));
    const double dt = value("--dt", ".001"), grain = value("--d50", ".00023");
    const double penalty = value("--penalty-time", ".01"), initialCs = value("--initial-cs", "0");
    const auto impact = get("--impact-model", "impulse");
    const double restitution = value("--restitution", ".5"), friction = value("--friction", ".3");
    if ((impact!="off" && impact!="impulse") || restitution<0 || restitution>1 || friction<0) {
        throw std::runtime_error("Use impact-model off or impulse, restitution in [0,1], and nonnegative friction");
    }
    const double grid = value("--refine", "1"), writeEvery = value("--write-interval", real(endTime).c_str());
    if (length<=0 || depth<=0 || burial<0 || endTime<dt || dt<=0 || grain<=0 || penalty<=0 ||
        initialCs<0 || initialCs>.1 || grid<1 || grid>3 || std::floor(grid)!=grid ||
        writeEvery<dt || writeEvery>endTime) throw std::runtime_error("Invalid physical or run parameters");
    const double segments = endTime / writeEvery;
    if (!std::isfinite(segments) || std::abs(segments-std::round(segments)) > 1e-10*segments) {
        throw std::runtime_error("Write interval must divide end time");
    }
    const auto inflow = readFlow(arguments.at("--flow"), endTime);
    const auto& flow = inflow.rows;
    auto faces = readHull(arguments.at("--stl"));
    const auto native = bounds(faces);
    const auto extent = subtract(native[1],native[0]);
    if (std::any_of(extent.begin(),extent.end(),[](double size){return size<=0;})) throw std::runtime_error("STL must have three positive extents");
    const double beam = value("--beam",real(length*extent[1]/extent[0]).c_str());
    const double height = value("--height",real(length*extent[2]/extent[0]).c_str());
    const double span = value("--domain-span",real(std::hypot(length,beam)).c_str());
    const double width = 1.7*span/(32*grid);
    const double thickness = value("--interface-thickness",real(width*.5).c_str());
    if (beam<=0 || height<=0 || span<=0 || burial>=height || height-burial>=depth || thickness<=0) throw std::runtime_error("Hull must intersect or touch the flat bed, with its top below the lid");
    const double a = std::remainder(yaw,360.)*std::acos(-1.)/180.;
    const double cosA = std::cos(a), sinA = std::sin(a);
    for (auto& face : faces) for (auto& p : face) {
        const double x = ((p[0]-native[0][0])/extent[0]-.5)*length;
        const double y = ((p[1]-native[0][1])/extent[1]-.5)*beam;
        p = {cx+cosA*x-sinA*y, cy+sinA*x+cosA*y, (p[2]-native[0][2])/extent[2]*height-burial};
    }
    verifyClosed(faces);
    const auto box = bounds(faces);
    if (std::max(std::abs(box[0][0]),std::abs(box[1][0]))+thickness>=.85*span ||
        std::max(std::abs(box[0][1]),std::abs(box[1][1]))+thickness>=.85*span ||
        box[1][2]+thickness>=depth) throw std::runtime_error("Hull and interface must fit within the refined region and below the lid");

    StagedCase stagedCase(output);
    const auto& stage = stagedCase.path();
    const auto u = [](const Flow& f){return tuple({f.speed,0,0});};
    const auto cs = [](const Flow& f){return real(f.concentration);};
    const auto k = [](const Flow& f){return real(1.5*std::pow(.05*f.speed,2));};
    const auto omega = [depth](const Flow& f){return real(std::sqrt(1.5)*.05*f.speed/(std::pow(.09,.25)*.07*depth));};
    const auto& first = flow.front();
    const std::string zeroGradient = "type zeroGradient;";
    write(stage/"system/blockMeshDict",mesh(span,depth,int(grid)));
    for (const char* file : {"fvSchemes","fvSolution","finite-area/faSchemes","finite-area/faSolution"}) {
        const auto dst = stage/"system"/file;
        fs::create_directories(dst.parent_path());
        fs::copy_file(templatePath/"system"/file,dst);
    }
    fs::create_directories(stage/"constant");
    for (const char* file : {"g","fvOptions","suspensionProperties","turbulenceProperties"}) {
        fs::copy_file(templatePath/"constant"/file,stage/"constant"/file);
    }
    write(stage/"system/finite-area/faMeshDefinition",dict("faMeshDefinition")+
          "polyMeshPatches (bed); boundary { inletFa {type patch; neighbourPolyPatch inlet;} outletFa {type patch; neighbourPolyPatch outlet;} sidesFa {type patch; neighbourPolyPatch sides;} } defaultPatch {name residual; type patch;}\n");
    write(stage/"system/controlDict",dict("controlDict")+
          "application shipwreckSedExnerFoam; startFrom startTime; startTime 0; stopAt endTime; endTime "+real(endTime)+
          "; deltaT "+real(dt)+"; maxDeltaT "+real(dt)+"; maxCo .2; adjustTimeStep true;\n"
          "writeControl adjustableRunTime; writeInterval "+real(writeEvery)+"; purgeWrite 2; writeFormat binary; writeCompression off; writePrecision 16; timePrecision 12; runTimeModifiable true;\n"
          "functions { bedShear {type wallShearStress; libs (fieldFunctionObjects); patches (bed); writeControl writeTime;} }\n");
    write(stage/"system/decomposeParDict",dict("decomposeParDict")+"numberOfSubdomains 6; method hierarchical; hierarchicalCoeffs {n (3 2 1); delta .001; order xyz;}\n");
    write(stage/"constant/shipwreckProperties",dict("shipwreckProperties")+"surfaceFile shipwreck.stl;\ncentre "+tuple({cx,cy,(box[0][2]+box[1][2])/2})+
          "; interfaceThickness "+real(thickness)+"; penaltyTime "+real(penalty)+";\nimpactEnabled "+
          (impact=="impulse"?"true":"false")+"; normalRestitution "+real(restitution)+"; frictionCoefficient "+real(friction)+";\n");
    write(stage/"constant/dynamicMeshDict",dict("dynamicMeshDict")+"dynamicFvMesh dynamicMotionSolverFvMesh; motionSolverLibs (\"libreferenceColumnMotion.so\"); solver referenceColumnMotion; referenceColumnMotionCoeffs {fixedTop "+real(depth)+";}\n");
    write(stage/"constant/transportProperties",dict("transportProperties")+"transportModel Newtonian; nu [0 2 -1 0 0 0 0] 1e-6; rhoF [1 -3 0 0 0 0 0] 1000; rhoS [1 -3 0 0 0 0 0] 2650; dS [0 1 0 0 0 0 0] "+real(grain)+"; CsMax [0 0 0 0 0 0 0] .57; Sc 1; reposeAngle 32;\n");
    write(stage/"constant/bedloadProperties",dict("bedloadProperties")+"sedimentBed on; bedMotion on; sedimentBedPatches (bed); writeSedimentBudget true; conservativeFixedPointProjection true; morphoAccFactor 1; ABorder 0; filterShields off; filterExner off; fixedFaPatches (inletFa outletFa); bedloadModel {type MeyerPeter; coefShields 1; avalanche on;} criticalShieldsModel {type Soulsby;}\n");
    volumeField(stage,"U","[0 1 -1 0 0 0 0]",u(first),table(flow,u),"type inletOutlet; inletValue uniform "+u(first)+"; value uniform "+u(first)+";","type movingWallVelocity; value uniform (0 0 0);",true);
    volumeField(stage,"p","[0 2 -2 0 0 0 0]","0",zeroGradient,"type fixedValue; value uniform 0;",zeroGradient);
    volumeField(stage,"k","[0 2 -2 0 0 0 0]",k(first),table(flow,k),zeroGradient,"type kqRWallFunction; value uniform "+k(first)+";");
    volumeField(stage,"omega","[0 0 -1 0 0 0 0]",omega(first),table(flow,omega),zeroGradient,"type omegaWallFunction; value uniform "+omega(first)+";");
    volumeField(stage,"nut","[0 2 -1 0 0 0 0]","0","type calculated; value uniform 0;","type calculated; value uniform 0;","type nutkWallFunction; value uniform 0;");
    volumeField(stage,"Cs","[0 0 0 0 0 0 0]",real(initialCs),table(flow,cs),"type inletOutlet; inletValue uniform 0; value uniform 0;",zeroGradient);
    volumeField(stage,"Ws","[0 1 -1 0 0 0 0]","(0 0 0)",zeroGradient,zeroGradient,zeroGradient,true);
    const std::string fixed = "type fixedValue; value uniform (0 0 0);";
    write(stage/"0/pointDisplacement",dict("pointDisplacement","pointVectorField")+"dimensions [0 1 0 0 0 0 0]; internalField uniform (0 0 0); boundaryField { inlet {"+fixed+"} outlet {"+fixed+"} bed {"+fixed+"} sides {type symmetry;} top {type symmetryPlane;} }\n");
    write(stage/"0/finite-area/qb",dict("qb","areaVectorField")+"dimensions [0 2 -1 0 0 0 0]; internalField uniform (0 0 0); boundaryField { inletFa {"+fixed+"} sidesFa {"+fixed+"} residual {"+fixed+"} outletFa {type zeroGradient;} }\n");
    write(stage/"inflow.csv", inflow.csv);
    write(stage/"shipwreck.foam","");
    write(stage/"scenario.json","{\n  \"length_m\": "+real(length)+", \"beam_m\": "+real(beam)+", \"height_m\": "+real(height)+
          ",\n  \"yaw_deg\": "+real(yaw)+", \"burial_m\": "+real(burial)+", \"centre_xy_m\": ["+real(cx)+", "+real(cy)+"],\n  \"depth_m\": "+real(depth)+
          ", \"domain_span_m\": "+real(span)+", \"end_time_s\": "+real(endTime)+", \"dt_s\": "+real(dt)+", \"cells\": "+std::to_string(66*52*28*int(grid*grid*grid))+
          ",\n  \"d50_m\": "+real(grain)+", \"interface_thickness_m\": "+real(thickness)+", \"penalty_time_s\": "+real(penalty)+
          ",\n  \"impact\": {\"model\": \""+impact+"\", \"normal_restitution\": "+real(restitution)+
          ", \"friction_coefficient\": "+real(friction)+"},\n  \"inlet_profile\": \"uniform in depth\", \"morphological_factor\": 1, \"experimental_validation_passed\": false\n}\n");
    writeHull(stage/"constant/triSurface/shipwreck.stl",faces);
    stagedCase.publish(output);
    std::cout << "Prepared " << output << ", bounds " << tuple(box[0]) << " to " << tuple(box[1]) << '\n';
}

int main(int argc,char** argv) {
    if (argc==2 && std::string(argv[1])=="--help") {
        std::cout << "shipwreck-case --case PATH --template PATH --stl BINARY_STL --flow CSV --length M --depth M --end S\n"
                     "  [--beam M --height M --yaw DEG --burial M --x M --y M --dt S --write-interval S]\n"
                     "  [--d50 M --initial-cs FRACTION --refine 1|2|3 --domain-span M --interface-thickness M --penalty-time S]\n"
                     "  [--impact-model impulse|off --restitution 0..1 --friction COEFFICIENT]\n";
        return 0;
    }
    try {
        const std::vector<std::string> keys{"--case","--template","--stl","--flow","--length","--depth","--end","--beam","--height","--yaw","--burial","--x","--y","--dt","--write-interval","--d50","--initial-cs","--refine","--domain-span","--interface-thickness","--penalty-time","--impact-model","--restitution","--friction"};
        std::map<std::string,std::string> arguments;
        for (int i=1; i<argc; i+=2) {
            if (i+1==argc || std::find(keys.begin(),keys.end(),argv[i])==keys.end() ||
                !arguments.emplace(argv[i],argv[i+1]).second) throw std::runtime_error("Unknown, duplicate or incomplete option; see --help");
        }
        for (std::size_t i=0; i<7; ++i) if (!arguments.count(keys[i])) throw std::runtime_error("Missing " + keys[i]);
        prepare(arguments);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
