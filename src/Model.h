#pragma once

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "GLSL.h"
#include "Program.h"

class Model {
public:
    Model();
    ~Model();

    // Core functionality
    bool loadModel(const std::string& path);
    void draw(std::shared_ptr<Program> shader);

    // Transform functions
    void setPosition(const glm::vec3& pos) { position = pos; }
    void setRotation(float angle, const glm::vec3& axis) { rotationAngle = angle; rotationAxis = axis; }
    void setScale(const glm::vec3& sc) { scale = sc; }

    // Animation functions
    void update(float deltaTime);
    void setAnimation(int animIndex);
    int getAnimationCount() const;
    std::string getAnimationName(int index) const;
    int getCurrentAnimation() const { return currentAnimation; }

    // Mesh information
    size_t getMeshCount() const { return meshes.size(); }
    bool hasTextures() const;

    // Material functions
    void setMeshMaterial(int meshIndex, const glm::vec3& ambient, const glm::vec3& diffuse,
        const glm::vec3& specular, float shininess);

private:
    // Data structures
    struct Vertex {
        glm::vec3 position;
        glm::vec3 normal;
        glm::vec2 texCoords;
        glm::vec3 tangent;
        glm::vec3 bitangent;
        glm::ivec4 boneIDs;
        glm::vec4 weights;
    };

    struct Texture {
        unsigned int id;
        std::string type;
        std::string path;
        bool embedded;
    };

    struct Material {
        glm::vec3 ambient;
        glm::vec3 diffuse;
        glm::vec3 specular;
        glm::vec3 emissive;
        float shininess;
        float opacity;
        bool hasTexture;
    };

    struct Mesh {
        std::vector<Vertex> vertices;
        std::vector<unsigned int> indices;
        std::vector<Texture> textures;
        Material material;
        GLuint VAO, VBO, EBO;
    };

    struct Bone {
        std::string name;
        int id;
        glm::mat4 offset;
    };

    struct Animation {
        std::string name;
        float duration;
        float ticksPerSecond;
        bool loop;
    };

    // Processing functions
    void processNode(aiNode* node, const aiScene* scene);
    Mesh processMesh(aiMesh* mesh, const aiScene* scene);
    std::vector<Texture> loadMaterialTextures(aiMaterial* mat, aiTextureType type, const std::string& typeName, const aiScene* scene);
    Material processMaterial(aiMaterial* material);
    void loadAnimations(const aiScene* scene);
    void updateBoneTransforms(float time, const aiScene* scene);

    // Texture loading utilities
    unsigned int loadTextureFromFile(const std::string& path);
    unsigned int loadEmbeddedTexture(const aiTexture* texture);

    // Member variables
    std::vector<Mesh> meshes;
    std::vector<Bone> bones;
    std::map<std::string, int> boneMap;
    std::vector<Animation> animations;
    std::vector<glm::mat4> boneTransforms;

    // Animation state
    int currentAnimation;
    float animationTime;

    // Transform state
    glm::vec3 position;
    float rotationAngle;
    glm::vec3 rotationAxis;
    glm::vec3 scale;

    // Asset management
    std::string directory;
    std::map<std::string, unsigned int> loadedTextures;

    // Assimp objects that must remain valid
    Assimp::Importer importer;
    const aiScene* scene;
};